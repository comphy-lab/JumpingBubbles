# Website Generation and Deployment

## Sitemap

```
.github/
├── README.md                  # This file
├── assets/                    # Website assets
│   ├── css/                   # CSS stylesheets
│   │   ├── academicons-1.7.0/ # Academic icons
│   │   ├── command-palette.css # Command palette styling
│   │   ├── custom_styles.css  # Custom styling
│   │   ├── fontello/          # Font icons
│   │   ├── fonts/             # Web fonts
│   │   ├── styles.css         # Main stylesheet
│   │   └── vendor.css         # Third-party styles
│   ├── js/                    # JavaScript files
│   │   ├── command-data.js    # Command palette data
│   │   ├── command-palette.js # Command palette functionality
│   │   ├── main.js            # Main JavaScript functionality
│   │   ├── search_db.json     # Search database
│   │   ├── shortcut-key.js    # Keyboard shortcuts
│   │   └── theme-toggle.js    # Theme switching
│   ├── custom_template.html   # HTML template for generated pages
│   ├── favicon/               # Favicon files
│   └── logos/                 # Logo images
├── scripts/                   # Build and deployment scripts
│   ├── build.sh               # Build script
│   ├── deploy.sh              # Deployment script
│   └── generate_docs.py       # Documentation generator
└── workflows/                 # GitHub Actions workflows
```

## Overview

This directory contains all the necessary files and scripts for generating and deploying the documentation of codes develped by CoMPhy Lab. The website is generated from source code files in the repository and is deployed at [comphy-lab.org/repositoryName](https://comphy-lab.org/repositoryName).

## Website Generation Process

The website generation process is handled by the `generate_docs.py` script, which converts source code files into HTML documentation. The process involves several steps:

1. **Source File Discovery**: The script scans the repository for source files in specified directories (`src-local`, `simulationCases`, `postProcess`).
2. **File Processing**: Each file is processed based on its type (C/C++, Python, Shell, Markdown).
3. **HTML Generation**: The processed content is converted to HTML using Pandoc with a custom template.
4. **Post-Processing**: The generated HTML is enhanced with additional features like code highlighting and navigation.
5. **Index Generation**: An index page is created from the README.md file.
6. **SEO Optimization**: Robots.txt and sitemap.xml files are generated for search engines.

### Key Components

#### 1. `generate_docs.py`

This Python script is the core of the documentation generation process. It performs the following functions:

- **Configuration Validation**: Ensures all required paths and files exist.
- **Source File Discovery**: Finds all source files in specified directories.
- **File Processing**: Processes different file types (C/C++, Python, Shell, Markdown) with specialized functions.
- **HTML Generation**: Uses Pandoc to convert processed content to HTML.
- **Post-Processing**: Enhances HTML with additional features.
- **Index Generation**: Creates an index page from README.md.
- **SEO Optimization**: Generates robots.txt and sitemap.xml.

#### 2. `build.sh`

This shell script runs the documentation generation process:

- Sets up the environment and paths.
- Runs the `generate_docs.py` script.
- Verifies that the documentation was generated successfully.

#### 3. `deploy.sh`

This shell script provides a local preview of the generated documentation:

- Verifies that the documentation exists.
- Starts a local web server on port 8000.
- Provides instructions for accessing the site.

#### 4. `custom_template.html`

This HTML template is used by Pandoc to generate the HTML pages:

- Defines the basic HTML structure.
- Includes meta tags for SEO.
- Links to CSS and JavaScript files.
- Provides placeholders for dynamic content.

## Assets

### CSS Files

- **styles.css**: Main stylesheet with comprehensive styling for the website.
- **custom_styles.css**: Custom styling specific to the CoMPhy Lab website.
- **command-palette.css**: Styling for the command palette feature.
- **vendor.css**: Third-party styles.

### JavaScript Files

- **main.js**: Main JavaScript functionality for the website.
- **command-palette.js**: Implements the command palette feature.
- **command-data.js**: Data for the command palette.
- **theme-toggle.js**: Handles theme switching (light/dark mode).
- **shortcut-key.js**: Implements keyboard shortcuts.
- **search_db.json**: Search database for the website.

## Deployment

The website is deployed at [comphy-lab.org](https://comphy-lab.org) using GitHub Pages. The deployment process involves:

1. Generating the documentation using `build.sh`.
2. Pushing the generated documentation to the `gh-pages` branch.
3. GitHub Pages serves the content from the `gh-pages` branch.

## Customization

To customize the website:

1. **Styling**: Modify the CSS files in `.github/assets/css/`.
2. **Functionality**: Modify the JavaScript files in `.github/assets/js/`.
3. **Template**: Modify the HTML template in `.github/assets/custom_template.html`.
4. **Generation Process**: Modify the Python script in `.github/scripts/generate_docs.py`.

## Troubleshooting

If the documentation generation fails:

1. Check the error messages from the `build.sh` script.
2. Verify that all required paths and files exist.
3. Ensure that Pandoc is installed and accessible.
4. Check that the Python script has the necessary permissions.

## Additional Information

- The website uses a custom template based on the Basilisk documentation.
- The documentation generation process is designed to handle C/C++, Python, Shell, and Markdown files.
- The website includes features like code highlighting, navigation, and search.
- SEO optimization is included to improve search engine visibility. 