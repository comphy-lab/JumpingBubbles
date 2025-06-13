#!/usr/bin/env python3
"""
Direct HTML Fix Script for Empty Anchors

This script directly removes empty anchor tags from HTML files that cause JavaScript syntax errors.
It uses regex-based string replacement without external dependencies, making it fast and lightweight.

Relationship to clean_html.py:
    - fix_empty_anchors.py: Fast, dependency-free approach using regular expressions. Works best
      with well-formed HTML where the empty anchor patterns are consistent. This script is ideal
      for quick cleanup or environments where installing dependencies is not possible.
    - clean_html.py: Uses BeautifulSoup for robust HTML parsing and sanitization. It can handle
      malformed HTML, complex nested structures, and content within script tags. If BeautifulSoup
      is not available, it falls back to the same regex approach used in this script.

Trade-offs:
    - Speed: fix_empty_anchors.py is generally faster as it avoids DOM parsing
    - Robustness: clean_html.py is more robust for complex or malformed HTML
    - Dependencies: fix_empty_anchors.py has no external dependencies
    - Script Tag Handling: clean_html.py better handles anchors in script tags

When to use which script:
    - Use fix_empty_anchors.py for: quick fixes, CI/CD environments with limited dependencies,
      well-formed HTML, performance-critical applications
    - Use clean_html.py for: complex HTML documents, when sanitization is important, or when
      working with HTML that may contain script tags with anchors

Dependencies:
    - None - relies only on standard library (re, os, glob)
    - Uses shared patterns from html_cleaning_patterns.py if available

Usage:
    python fix_empty_anchors.py [options] <path>
"""

import sys
import re
import os
import glob
import argparse

# Try to import shared patterns, fall back to local patterns if not available
try:
    from html_cleaning_patterns import EMPTY_ANCHOR_PATTERNS, apply_empty_anchor_cleanup
    USES_SHARED_PATTERNS = True
except ImportError:
    # If shared patterns aren't available, define patterns locally
    USES_SHARED_PATTERNS = False
    # Define patterns here for backwards compatibility
    EMPTY_ANCHOR_PATTERNS = [
        # Combined pattern for id and href attributes in either order with optional newlines
        r'(?i)<a\s+id=[\'"]?([^\s>]*)[\'"]?\s+href=[\'"]?#[\'"]?\s*>\s*(?:\n\s*)*</a>',
        r'(?i)<a\s+href=[\'"]?#[\'"]?\s+id=[\'"]?([^\s>]*)[\'"]?\s*>\s*(?:\n\s*)*</a>',
        
        # Single pattern for id attribute only with optional newlines
        r'(?i)<a\s+id=[\'"]?([^\s>]*)[\'"]?\s*>\s*(?:\n\s*)*</a>',
        
        # Single pattern for href='#' only with optional newlines
        r'(?i)<a\s+href=[\'"]?#[\'"]?\s*>\s*(?:\n\s*)*</a>',
        
        # Combined pattern for unquoted attributes in either order
        r'(?i)<a\s+(?:id=([^\s>]*)\s+href=#|href=#\s+id=([^\s>]*))\s*>\s*(?:\n\s*)*</a>'
    ]
    
    # Implement apply_empty_anchor_cleanup locally for consistency
    def apply_empty_anchor_cleanup(content):
        """
        Apply all empty anchor cleanup patterns to the content.
        Local implementation for when the shared module is not available.
        
        Args:
            content: HTML content to clean
            
        Returns:
            str: Cleaned HTML content with empty anchors removed
        """
        result = content
        for pattern in EMPTY_ANCHOR_PATTERNS:
            result = re.sub(pattern, '', result)
        return result

def fix_html_file(file_path, verbose=False, dry_run=False):
    """
    Removes empty anchor tags from an HTML file using direct string replacement.

    Args:
        file_path: Path to the HTML file to clean
        verbose: If True, prints detailed information
        dry_run: If True, shows changes without writing to file

    Returns:
        int: Number of replacements made
    """
    try:
        # Read the file
        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read()

        # Store original content length for comparison
        initial_content_length = len(content)

        if verbose:
            print(f"Processing file: {file_path}")
        
        # Apply anchor cleanup using shared patterns if available,
        # or fall back to local implementation
        if USES_SHARED_PATTERNS and 'apply_empty_anchor_cleanup' in globals():
            # Use the imported shared function
            modified_content = apply_empty_anchor_cleanup(content)
        else:
            # Apply all patterns directly
            modified_content = content
            for pattern in EMPTY_ANCHOR_PATTERNS:
                modified_content = re.sub(pattern, '', modified_content)

        # Calculate approximate number of replacements
        final_content_length = len(modified_content)
        chars_removed = initial_content_length - final_content_length
        replacements = chars_removed // 20  # Approximate size of each anchor tag

        # Only proceed if changes were made
        if replacements > 0:
            if verbose:
                print(f"  Found approximately {replacements} empty anchor tags")

            if dry_run:
                print(f"[DRY RUN] Would fix {file_path}: {replacements} empty anchor tags")
            else:
                with open(file_path, 'w', encoding='utf-8') as f:
                    f.write(modified_content)
                print(f"Fixed {file_path}: removed approximately {replacements} empty anchor tags")
        elif verbose:
            print(f"  No empty anchors found in {file_path}")
        
        return replacements
        
    except Exception as e:
        print(f"Error fixing {file_path}: {e}")
        return 0

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description="Fix empty anchor tags in HTML files that cause JavaScript syntax errors.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "path",
        help="Path to an HTML file or directory containing HTML files"
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose output"
    )
    parser.add_argument(
        "-d", "--dry-run",
        action="store_true",
        help="Show what would be fixed without making changes"
    )

    args = parser.parse_args()

    path = args.path
    verbose = args.verbose
    dry_run = args.dry_run

    if dry_run:
        print("Running in dry-run mode - no changes will be made")

    if os.path.isfile(path):
        # Fix a single file
        fix_html_file(path, verbose=verbose, dry_run=dry_run)
    elif os.path.isdir(path):
        # Fix all HTML files in the directory and subdirectories
        if verbose:
            print(f"Searching for HTML files in {path} and subdirectories...")

        html_files = glob.glob(os.path.join(path, '**', '*.html'), recursive=True)
        total_files = len(html_files)

        if verbose:
            print(f"Found {total_files} HTML files")

        fixed_files = 0
        total_replacements = 0

        for file in html_files:
            replacements = fix_html_file(file, verbose=verbose, dry_run=dry_run)
            if replacements > 0:
                fixed_files += 1
                total_replacements += replacements

        action = "Would fix" if dry_run else "Fixed"
        print(f"\nSummary: {action} {fixed_files} out of {total_files} files, removing approximately {total_replacements} empty anchor tags")
    else:
        print(f"Error: {path} is not a valid file or directory")
        sys.exit(1)

if __name__ == '__main__':
    main()