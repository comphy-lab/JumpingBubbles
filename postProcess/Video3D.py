#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Basilisk C Snapshot Visualization Script

This script processes Basilisk C simulation output files to generate 3D visualizations
of multiphase flow simulations. It reads cell and facet data from intermediate 
snapshot files, creates 3D meshes using PyVista, and saves rendered images for
animation production.

The script uses multiprocessing to handle large datasets efficiently and includes
memory management strategies to prevent out-of-memory errors during batch processing.

Usage:
    python visualize_snapshots.py

Dependencies:
    - numpy: Array processing
    - pyvista: 3D visualization and mesh handling
    - multiprocessing: Parallel processing
    - Basilisk C executables: getCells_bottomPlate, getFacets3D
"""

import subprocess as sp
import numpy as np
import pyvista as pv
import os
import multiprocessing
import gc  # Garbage collection

# ===============================
# Configuration and Settings
# ===============================

# Camera configuration for consistent visualization
CAMERA_CONFIG = {
    "position": (7.75, 3.50, -2.50),
    "focal_point": (0.0, 0.0, 0.0),
    "up_vector": (-0.30, 0.95, 0.025),
    "parallel_scale": 2.2
}

# Reflection planes for symmetry operations
REFLECTION_NORMALS = {
    'x': [1, 0, 0],
    'z': [0, 0, 1]
}

# ===============================
# Utility Functions
# ===============================

def parse_vertex(vertex_str):
    """
    Parse a vertex string into a numpy array.

    Converts space-separated string coordinates into a numpy array
    for 3D point representation.

    Args:
        vertex_str (str): Space-separated vertex coordinates (e.g., "1.0 2.0 3.0")

    Returns:
        np.ndarray: 3D coordinate array

    Example:
        >>> parse_vertex("1.0 2.0 3.0")
        array([1., 2., 3.])
    """
    return np.fromstring(vertex_str, sep=' ')


def run_process(command):
    """
    Execute a Basilisk C utility and capture its output.

    Runs external Basilisk processing tools and returns the stderr output,
    which contains the processed geometric data.

    Args:
        command (list): Command and arguments to execute

    Returns:
        str: Decoded stderr output containing geometric data

    Raises:
        subprocess.CalledProcessError: If the external command fails

    Note:
        Basilisk C utilities output data to stderr by convention
    """
    p = sp.Popen(command, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    if p.returncode != 0:
        error_msg = stdout.decode("utf-8") if stdout else stderr.decode("utf-8")
        raise sp.CalledProcessError(p.returncode, command[0], output=error_msg)
    return stderr.decode("utf-8").strip()


# ===============================
# Mesh Creation Functions
# ===============================

def create_polydata(points, lines_array, faces_array):
    """
    Create a PyVista PolyData object from points and connectivity arrays.

    Constructs a complete mesh representation with vertices, edges, and faces
    for visualization purposes.

    Args:
        points (list): List of 3D points defining vertices
        lines_array (np.ndarray): Connectivity array for line elements
        faces_array (np.ndarray): Connectivity array for face elements

    Returns:
        pv.PolyData: Complete mesh object ready for visualization
    """
    polydata = pv.PolyData()
    polydata.points = np.array(points)
    polydata.lines = lines_array
    polydata.faces = faces_array
    return polydata


def process_cells(cell_data):
    """
    Process cell data from Basilisk output into PyVista meshes.

    Parses the cell structure data and creates individual cell meshes
    representing the computational grid elements.

    Args:
        cell_data (str): Raw cell data from getCells_bottomPlate output

    Returns:
        list: List of PyVista PolyData objects representing cells

    Note:
        Assumes rectangular cells with standard connectivity pattern
    """
    # Standard connectivity for rectangular cells
    lines_array = np.array([[2, 0, 1], [2, 1, 2], [2, 2, 3], [2, 3, 0]], dtype=np.int32)
    faces_array = np.array([[4, 0, 1, 2, 3]], dtype=np.int32)
    
    # Parse cell data blocks separated by double newlines
    cells = [[parse_vertex(point) for point in cell.split('\n')] 
             for cell in cell_data.split('\n\n')]
    
    return [create_polydata(cell, lines_array, faces_array) for cell in cells]


def process_facets(facet_data):
    """
    Process facet data into a single unified mesh.

    Converts facet vertex data into an indexed mesh representation,
    eliminating duplicate vertices for efficiency. This is typically
    used for interface reconstruction in multiphase flows.

    Args:
        facet_data (str): Raw facet data from getFacets3D output

    Returns:
        pv.PolyData: Unified mesh containing all facets

    Performance:
        Uses dictionary-based vertex indexing to avoid duplicate points,
        reducing memory usage for large datasets.
    """
    vertex_to_index = {}
    points, cells = [], []
    
    for facet in facet_data.split('\n\n'):
        cell = []
        for vertex_str in facet.split('\n'):
            vertex = parse_vertex(vertex_str)
            vertex_tuple = tuple(vertex)  # Convert to hashable type
            
            # Index unique vertices
            if vertex_tuple not in vertex_to_index:
                vertex_to_index[vertex_tuple] = len(points)
                points.append(vertex)
            cell.append(vertex_to_index[vertex_tuple])
        
        # Build VTK-style connectivity array
        cells.extend([len(cell)] + cell)
    
    return pv.PolyData(np.array(points, dtype=np.float64), 
                      faces=np.array(cells, dtype=np.int32))


def reflect_mesh(mesh, normals):
    """
    Apply reflection symmetry operations to a mesh.

    Creates reflected copies of the mesh across specified planes to
    reconstruct the full domain from symmetry-reduced simulations.

    Args:
        mesh (pv.PolyData): Original mesh to reflect
        normals (dict): Dictionary mapping labels to normal vectors

    Returns:
        pv.PolyData: Combined mesh including all reflections

    Note:
        Reflections are applied sequentially and merged with the original
    """
    for normal in normals.values():
        mesh = mesh.merge(mesh.reflect(normal))
    return mesh


# ===============================
# Main Processing Functions
# ===============================

def get_grid(filename):
    """
    Extract and process computational grid from Basilisk snapshot.

    Reads cell data from a snapshot file and constructs the visualization
    mesh for the computational grid, applying symmetry operations.

    Args:
        filename (str): Path to Basilisk snapshot file

    Returns:
        pv.PolyData: Complete grid mesh with reflections applied

    Note:
        Requires getCells_bottomPlate executable in the working directory
    """
    cell_data = run_process(["./getCells_bottomPlate", filename])
    cells = process_cells(cell_data)
    mesh = pv.MultiBlock(cells).combine()
    return reflect_mesh(mesh, REFLECTION_NORMALS)


def get_facets_3d(filename):
    """
    Extract and process 3D interface facets from Basilisk snapshot.

    Reads facet data representing fluid interfaces and constructs
    the visualization mesh with symmetry operations applied.

    Args:
        filename (str): Path to Basilisk snapshot file

    Returns:
        pv.PolyData: Complete interface mesh with reflections

    Note:
        Requires getFacets3D executable in the working directory
    """
    facet_data = run_process(["./getFacets3D", filename])
    mesh = process_facets(facet_data)
    return reflect_mesh(mesh, REFLECTION_NORMALS)


def process_and_save_image(t, base_filename, image_folder):
    """
    Process a single time step and generate visualization image.

    Reads simulation data for a specific time, creates 3D visualization
    with grid and interface meshes, and saves the rendered image.

    Args:
        t (float): Simulation time
        base_filename (str): Template filename with format specifier
        image_folder (str): Output directory for images

    Performance:
        Includes explicit garbage collection to manage memory usage
        during batch processing of large datasets.

    Note:
        Skips processing if output image already exists or input file
        is missing. Uses off-screen rendering for headless operation.
    """
    filename = base_filename % t
    image_filename = f"{image_folder}/{int(t*1e4):06d}.png"
    
    # Check file existence
    if not os.path.exists(filename):
        print(f"File {filename} does not exist")
        return
    if os.path.exists(image_filename):
        print(f"Image {image_filename} already exists")
        return
    
    print(f"Processing t = {t}")

    # Process mesh data
    cells = get_grid(filename)
    poly_data = get_facets_3d(filename)

    # Configure off-screen renderer
    plotter = pv.Plotter(off_screen=True)

    try:
        # Add mesh components with different styles
        plotter.add_mesh(cells, color='grey')  # Solid grid
        plotter.add_mesh(cells, color='black', style='wireframe')  # Grid edges
        plotter.add_mesh(poly_data, color='orange')  # Interface

        # Apply camera configuration
        plotter.camera.position = CAMERA_CONFIG["position"]
        plotter.camera.focal_point = CAMERA_CONFIG["focal_point"]
        plotter.camera.up = CAMERA_CONFIG["up_vector"]
        plotter.camera.parallel_scale = CAMERA_CONFIG["parallel_scale"]

        # Add timestamp annotation
        plotter.add_text(f"t = {t:.2f}", position='upper_right', 
                        font_size=15, color='black', font_file='cmunrm.ttf')
        
        # Render and save
        plotter.screenshot(image_filename)
        plotter.close()
        
    except Exception as e:
        print(f"Error processing t = {t}: {str(e)}")
        
    finally:
        # Explicit memory cleanup for large datasets
        del cells, poly_data, plotter
        gc.collect()


def process_single_time_step(args):
    """
    Wrapper function for multiprocessing pool.

    Unpacks arguments and calls the main processing function,
    ensuring proper garbage collection after each process.

    Args:
        args (tuple): Packed arguments for process_and_save_image

    Note:
        Additional garbage collection helps prevent memory accumulation
        in long-running batch processes.
    """
    process_and_save_image(*args)
    gc.collect()


# ===============================
# Main Execution
# ===============================

def main():
    """
    Main function that drives the script execution.

    Sets up the processing pipeline for batch visualization of Basilisk
    snapshots, using multiprocessing to handle multiple time steps in
    parallel. Creates output directory and manages the processing pool.
    """
    # Configuration
    base_filename = "intermediate/snapshot-%5.4f"
    image_folder = "Video"
    
    # Time stepping parameters
    time_step = 0.01
    end_time = 8
    
    # Create output directory
    if not os.path.exists(image_folder):
        os.makedirs(image_folder)

    # Generate time steps
    time_steps = np.arange(0, end_time + time_step, time_step)
    args = [(t, base_filename, image_folder) for t in time_steps]

    # Process in parallel with limited pool size for memory management
    # TODO: Make number of processes configurable via command line
    with multiprocessing.Pool(processes=2) as pool:
        pool.map(process_single_time_step, args)


if __name__ == '__main__':
    main()