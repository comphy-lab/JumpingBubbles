/**
 * Command Palette functionality for CoMPhy Lab website
 * This file contains all the functionality for the command palette
 */

// Use centralized DEBUG flag from config
const DEBUG = window.CONFIG ? window.CONFIG.DEBUG : false;

// Make the command palette opening function globally available
window.openCommandPalette = function() {
  const palette = document.getElementById('simple-command-palette');
  if (palette) {
    palette.style.display = 'block';
    const input = document.getElementById('command-palette-input');
    if (input) {
      input.value = '';
      input.focus();
      if (typeof renderCommandResults === 'function') {
        renderCommandResults('');
      }
    }
  }
};

/**
 * Renders the command palette results filtered by the provided query.
 *
 * Filters available commands by title or section, groups them by section, and displays them in the results container. If the query is at least three characters and an asynchronous search function is available, performs an additional database search and updates the results with any matches. Displays a "No commands found" message if no results are available.
 *
 * @param {string} query - The search string used to filter and search commands.
 */
function renderCommandResults(query) {
  const resultsContainer = document.getElementById('command-palette-results');
  if (!resultsContainer) return;
  
  // Clear results
  resultsContainer.innerHTML = '';
  
  // Get commands
  const commands = window.commandData || [];
  
  // Filter commands based on query
  const filteredCommands = query ? 
    commands.filter(cmd => 
      cmd.title.toLowerCase().includes(query.toLowerCase()) ||
      (cmd.section && cmd.section.toLowerCase().includes(query.toLowerCase()))
    ) : 
    commands;
  
  // Group by section
  const sections = {};
  filteredCommands.forEach(cmd => {
    if (!sections[cmd.section]) {
      sections[cmd.section] = [];
    }
    sections[cmd.section].push(cmd);
  });
  
  // If query is at least 3 characters, search the database as well
  if (query && query.length >= 3 && typeof window.searchDatabaseForCommandPalette === 'function') {
    // Capture the current query to avoid stale updates
    const capturedQuery = query;
    
    // We'll use a promise to handle the async search
    window.searchDatabaseForCommandPalette(capturedQuery).then(searchResults => {
      // Get the current input value
      const currentInputValue = document.getElementById('command-palette-input')?.value || '';
      
      // Only update UI if the captured query matches the current input value
      if (capturedQuery === currentInputValue && searchResults && searchResults.length > 0) {
        // Add search results to sections
        sections['Search Results'] = searchResults;
        
        // Re-render the UI with search results
        renderSections(sections, resultsContainer);
      }
    }).catch(err => {
      if (DEBUG) {
        console.error('Error searching database:', err);
      }
    });
  }
  
  // Render the sections we have now (this will be called immediately, and again if search results come in)
  renderSections(sections, resultsContainer);
  
  // Show message if no results
  if (Object.keys(sections).length === 0) {
    const noResults = document.createElement('div');
    noResults.className = 'command-palette-no-results';
    noResults.textContent = 'No commands found';
    resultsContainer.appendChild(noResults);
  }
}

/**
 * Renders grouped command sections into a container element for the command palette UI.
 *
 * Each section is displayed with its title and a list of commands, including icons, titles, and optional excerpts. Command icons that use inline HTML (such as Font-Awesome) are rendered as raw HTML if their source is constant; otherwise, icons are rendered as plain text. Clicking a command triggers its handler and hides the palette.
 *
 * @param {Object} sections - An object mapping section titles to arrays of command objects.
 * @param {HTMLElement} container - The DOM element where the sections will be rendered.
 */
function renderSections(sections, container) {
  // Clear container first
  container.innerHTML = '';
  
  // Create DOM elements for results
  Object.keys(sections).forEach(section => {
    const sectionEl = document.createElement('div');
    sectionEl.className = 'command-palette-section';
    
    // Add additional class for search results section to style it properly
    if (section === 'Search Results') {
      sectionEl.classList.add('command-palette-search-results-section');
    }
    
    const sectionTitle = document.createElement('div');
    sectionTitle.className = 'command-palette-section-title';
    sectionTitle.textContent = section;
    sectionEl.appendChild(sectionTitle);
    
    const commandsList = document.createElement('div');
    commandsList.className = 'command-palette-commands';
    
    // Simple HTML sanitizer that only allows specific tags and attributes
    function sanitizeHTML(html) {
      // Only allow these tags and attributes
      const allowedTags = ['i', 'span', 'svg', 'path'];
      const allowedAttrs = ['class', 'viewBox', 'd', 'fill', 'stroke', 'stroke-width', 'xmlns'];
      
      // Create a temporary element to parse the HTML
      const temp = document.createElement('div');
      temp.innerHTML = html;
      
      // Function to sanitize a node and its children
      function sanitizeNode(node) {
        // For element nodes
        if (node.nodeType === 1) {
          // If it's not an allowed tag, replace it with its text content
          if (!allowedTags.includes(node.tagName.toLowerCase())) {
            return document.createTextNode(node.textContent);
          }
          
          // Remove any attributes that aren't allowed
          Array.from(node.attributes).forEach(attr => {
            if (!allowedAttrs.includes(attr.name.toLowerCase())) {
              node.removeAttribute(attr.name);
            }
          });
          
          // Sanitize all child nodes
          Array.from(node.childNodes).forEach(child => {
            const sanitizedChild = sanitizeNode(child);
            if (sanitizedChild !== child) {
              node.replaceChild(sanitizedChild, child);
            }
          });
        }
        
        return node;
      }
      
      // Sanitize the temporary element and all its children
      sanitizeNode(temp);
      
      // Return the sanitized HTML
      return temp.innerHTML;
    }

    sections[section].forEach(cmd => {
      const cmdEl = document.createElement('div');
      cmdEl.className = 'command-palette-command';
      
      // Create icon element separately with sanitized HTML
      const iconEl = document.createElement('div');
      iconEl.className = 'command-palette-icon';
      
      // Render Font-Awesome & other inline-HTML icons -> sanitize to prevent XSS risks
      if (cmd.icon && typeof cmd.icon === 'string') {
        // Strict validation: only allow Font Awesome icons that match the expected pattern
        const iconStr = cmd.icon.trim();
        // Only allow Font Awesome icons with specific classes or SVG paths
        if (iconStr.startsWith('<i') && 
            (/class=\"(fa-|ai )/.test(iconStr) || iconStr.includes('fa-solid') ||
             iconStr.includes('fa-brands') || iconStr.includes('fa-regular') ||
             iconStr.includes('ai ai-'))) {
          // Adding extra validation for classes
          const validIconPattern = /^<i class=\"(fa|fas|far|fal|fab|ai)[\s-][\w\d\s-]+"><\/i>$/;
          if (validIconPattern.test(iconStr)) {
            iconEl.innerHTML = sanitizeHTML(iconStr);
          } else {
            // Fallback to displaying as text if pattern doesn't match exactly
            iconEl.textContent = '🔍'; // Default search icon as fallback
            if (DEBUG) {
              console.warn('Invalid icon format detected:', iconStr);
            }
          }
        } else {
          // Just use text for non-icon strings
          iconEl.textContent = iconStr;
        }
      } else {
        // Default icon if none provided
        iconEl.textContent = '';
      }
      
      // Create title element
      const titleEl = document.createElement('div');
      titleEl.className = 'command-palette-title';
      // Use our sanitizer to prevent XSS if available, otherwise set textContent directly
      if (window.htmlSanitizer && typeof window.htmlSanitizer.setTextContent === 'function') {
        window.htmlSanitizer.setTextContent(titleEl, cmd.title);
      } else {
        // Fallback to direct textContent setting if sanitizer is not available
        titleEl.textContent = cmd.title || '';
      }

      // Build command element
      cmdEl.appendChild(iconEl);
      cmdEl.appendChild(titleEl);

      // Add excerpt for search results if available
      if (cmd.excerpt) {
        const excerptEl = document.createElement('div');
        excerptEl.className = 'command-palette-excerpt';
        // Use our sanitizer to prevent XSS
        const excerptText = cmd.excerpt.substring(0, 120) + (cmd.excerpt.length > 120 ? '...' : '');
        
        // Use our sanitizer to prevent XSS if available, otherwise set textContent directly
        if (window.htmlSanitizer && typeof window.htmlSanitizer.setTextContent === 'function') {
          window.htmlSanitizer.setTextContent(excerptEl, excerptText);
        } else {
          // Fallback to direct textContent setting if sanitizer is not available
          excerptEl.textContent = excerptText || '';
        }
        
        cmdEl.appendChild(excerptEl);
      }
      
      cmdEl.addEventListener('click', function(e) {
        if (typeof cmd.handler === 'function') {
          document.getElementById('simple-command-palette').style.display = 'none';
          cmd.handler();
        }
      });
      
      commandsList.appendChild(cmdEl);
    });
    
    sectionEl.appendChild(commandsList);
    container.appendChild(sectionEl);
  });
}

/**
 * Initializes the command palette UI and search functionality on page load.
 *
 * Sets up event listeners for keyboard shortcuts, input handling, and UI interactions. Prefetches the search database and initializes fuzzy search if available. Enables keyboard navigation, command execution, and closing the palette via backdrop or Escape key.
 *
 * @remark If the search database cannot be fetched, search functionality will be unavailable, but the command palette UI will still operate with local commands.
 */
function initCommandPalette() {
  // Ensure search database is preloaded for command palette search functionality
  // Try to prefetch the search database if it exists
  // Get base URL from meta tag to support GitHub Pages subfolders
  const baseUrlMeta = document.querySelector('meta[name="base-url"]');
  const baseUrl = baseUrlMeta ? baseUrlMeta.getAttribute('content') : '';
  const searchDbUrl = baseUrl ? `${baseUrl}/assets/js/search_db.json` : '/assets/js/search_db.json';
  
  fetch(searchDbUrl).then(response => {
    if (response.ok) {
      return response.json();
    }
    throw new Error('Search database not found');
  }).then(data => {
    if (DEBUG) {
      console.log('Search database prefetched for command palette');
    }
    window.searchData = data;
    
    // Initialize Fuse.js with weighted keys
    if (typeof Fuse === 'function') {
      window.searchFuse = new Fuse(data, {
        keys: [
          { name: 'title', weight: 0.7 },
          { name: 'content', weight: 0.2 },
          { name: 'tags', weight: 0.1 },
          { name: 'categories', weight: 0.1 }
        ],
        includeScore: true,
        threshold: 0.4
      });
    } else {
      console.warn('Fuse.js not found – command-palette search disabled');
    }
  }).catch(err => {
    if (DEBUG) {
      console.warn('Could not prefetch search database for command palette:', err.message);
    }
  });

  // Set up backdrop click to close
  const backdrop = document.querySelector('.simple-command-palette-backdrop');
  if (backdrop) {
    backdrop.addEventListener('click', function() {
      document.getElementById('simple-command-palette').style.display = 'none';
    });
  }
  
  // Set up input handler
  const input = document.getElementById('command-palette-input');
  if (input) {
    input.addEventListener('input', function() {
      renderCommandResults(this.value);
    });
    
    input.addEventListener('keydown', function(e) {
      if (e.key === 'Escape') {
        document.getElementById('simple-command-palette').style.display = 'none';
      } else if (e.key === 'Enter') {
        const selectedCommand = document.querySelector('.command-palette-command.selected');
        if (selectedCommand) {
          selectedCommand.click();
        } else {
          const firstCommand = document.querySelector('.command-palette-command');
          if (firstCommand) {
            firstCommand.click();
          }
        }
      } else if (e.key === 'ArrowDown' || e.key === 'ArrowUp') {
        e.preventDefault();
        
        const commands = Array.from(document.querySelectorAll('.command-palette-command'));
        if (commands.length === 0) return;
        
        const currentSelected = document.querySelector('.command-palette-command.selected');
        let nextIndex = 0;
        
        if (currentSelected) {
          const currentIndex = commands.indexOf(currentSelected);
          currentSelected.classList.remove('selected');
          
          if (e.key === 'ArrowDown') {
            nextIndex = (currentIndex + 1) % commands.length;
          } else {
            nextIndex = (currentIndex - 1 + commands.length) % commands.length;
          }
        } else {
          nextIndex = e.key === 'ArrowDown' ? 0 : commands.length - 1;
        }
        
        commands[nextIndex].classList.add('selected');
        
        // Ensure the selected element is visible in the scroll view
        commands[nextIndex].scrollIntoView({
          behavior: 'smooth',
          block: 'nearest'
        });
      }
    });
  }
  
  // Register command palette keyboard shortcut
  document.addEventListener('keydown', function(e) {
    if ((e.metaKey || e.ctrlKey) && e.key === 'k') {
      e.preventDefault();
      window.openCommandPalette();
    }
  });
  
  // Make command palette button work
  const commandPaletteBtn = document.getElementById('command-palette-btn');
  if (commandPaletteBtn) {
    commandPaletteBtn.addEventListener('click', function(e) {
      e.preventDefault();
      window.openCommandPalette();
    });
  }
}

// Run initialization when DOM is loaded
document.addEventListener('DOMContentLoaded', function() {
  // Initialize the command palette
  initCommandPalette();
  
  // Show appropriate shortcut text based on platform
  // Use modern platform detection methods
  let isMac = false;
  if (navigator.userAgentData && navigator.userAgentData.platform) {
    // Use the newer userAgentData API if available
    isMac = navigator.userAgentData.platform.includes('macOS');
  } else {
    // Fall back to userAgent string for older browsers
    isMac = /(Mac|iPhone|iPod|iPad)/i.test(navigator.userAgent);
  }
  
  document.querySelectorAll('.mac-theme-text').forEach(el => {
    el.style.display = isMac ? 'inline' : 'none';
  });
  document.querySelectorAll('.default-theme-text').forEach(el => {
    el.style.display = isMac ? 'none' : 'inline';
  });
  
  // Set the appropriate shortcut hint based on platform
  const shortcutHint = document.getElementById('command-palette-shortcut');
  if (shortcutHint) {
    shortcutHint.textContent = isMac ? '⌘K' : 'Ctrl+K';
  }
  
  // Ensure command palette button works correctly
  const commandPaletteBtn = document.getElementById('command-palette-btn');
  if (commandPaletteBtn) {
    if (DEBUG) {
      console.log('Command palette button initialized with new styling');
    }
    
    // Make sure the button retains focus styles
    commandPaletteBtn.addEventListener('focus', function() {
      this.classList.add('focused');
    });
    
    commandPaletteBtn.addEventListener('blur', function() {
      this.classList.remove('focused');
    });
  }
});

// Make functions available globally
window.renderCommandResults = renderCommandResults;
window.renderSections = renderSections;

// Use the shared search database function from search-helper.js
window.searchDatabaseForCommandPalette = async function(query) {
  if (window.searchHelper && typeof window.searchHelper.searchDatabaseForCommandPalette === 'function') {
    return window.searchHelper.searchDatabaseForCommandPalette(query);
  }
  return [];
};
