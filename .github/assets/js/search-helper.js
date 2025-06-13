/**
 * Shared search helper functions for website command palette
 * This file contains the shared search functionality used by both command-data.js and command-palette.js
 */

// Initialize the window.searchHelper namespace if it doesn't exist
window.searchHelper = window.searchHelper || {};

// Store the promise of search database and Fuse creation for memoization
window.searchHelper.fuseInitPromise = null;

// Initialize the search database and Fuse instance (memoized)
window.searchHelper.initializeSearchFuse = async function() {
  // If the initialization is already in progress, return the existing promise
  if (window.searchHelper.fuseInitPromise) {
    return window.searchHelper.fuseInitPromise;
  }

  // If Fuse is already initialized, return a resolved promise
  if (window.searchFuse) {
    return Promise.resolve();
  }

  // Create a new initialization promise
  window.searchHelper.fuseInitPromise = (async () => {
    try {
      // If we already have search data but no Fuse object
      if (window.searchData) {
        // Check if Fuse is available before creating a new instance
        if (typeof Fuse !== 'undefined') {
          window.searchFuse = new Fuse(window.searchData, {
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
          console.warn('Fuse.js is not loaded. Search functionality will be unavailable.');
          // Create a flag to indicate Fuse is not available
          window.searchFuseUnavailable = true;
        }
        return;
      }

      // Get base URL from meta tag to support GitHub Pages subfolders
      const baseUrlMeta = document.querySelector('meta[name="base-url"]');
      const baseUrl = baseUrlMeta ? baseUrlMeta.getAttribute('content') : '';
      const searchDbUrl = baseUrl ? `${baseUrl}/assets/js/search_db.json` : '/assets/js/search_db.json';

      const response = await fetch(searchDbUrl);
      if (!response.ok) {
        console.warn(`No search database found (${response.status})`);
        throw new Error(`Failed to fetch search database: ${response.status}`);
      }

      const searchData = await response.json();
      if (!searchData || !Array.isArray(searchData)) {
        console.warn('Search database has invalid format');
        throw new Error('Search database has invalid format');
      }

      window.searchData = searchData;

      // Check if Fuse is available before creating a new instance
      if (typeof Fuse !== 'undefined') {
        window.searchFuse = new Fuse(searchData, {
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
        console.warn('Fuse.js is not loaded. Search functionality will be unavailable.');
        // Create a flag to indicate Fuse is not available
        window.searchFuseUnavailable = true;
      }
    } catch (e) {
      console.error('Error initializing search:', e);
      throw e;
    }
  })();

  return window.searchHelper.fuseInitPromise;
};

// Shared search database function for command palette
window.searchHelper.searchDatabaseForCommandPalette = async function(query) {
  // Only perform search if query is at least 3 characters long
  if (!query || query.length < 3) {
    return [];
  }

  console.log('Searching database for:', query);

  try {
    // Initialize Fuse if not already initialized
    await window.searchHelper.initializeSearchFuse();
    
    // Check if Fuse is unavailable after initialization attempt
    if (window.searchFuseUnavailable) {
      console.warn('Search attempted but Fuse.js is not available.');
      return [];
    }

    // Perform the search with the initialized Fuse instance
    if (window.searchFuse) {
      try {
        const results = window.searchFuse.search(query);

        // Sort results by priority first, then by Fuse.js score
        // Lower priority number = higher priority (1 is highest, 5 is lowest)
        const sortedResults = results.sort((a, b) => {
          // First compare by priority
          const priorityA = a.item.priority || a.item.priority === 0 ? a.item.priority : 5; // Default to lowest priority if not specified
          const priorityB = b.item.priority || b.item.priority === 0 ? b.item.priority : 5;

          if (priorityA !== priorityB) {
            return priorityA - priorityB; // Lower priority number comes first
          }

          // If priorities are equal, use Fuse.js score (lower score = better match)
          return a.score - b.score;
        });

        // Helper function for basic HTML escaping if htmlSanitizer is not available
        const escapeHtml = (str) => {
          if (!str) return '';
          return String(str)
            .replace(/&/g, '&amp;')
            .replace(/</g, '&lt;')
            .replace(/>/g, '&gt;')
            .replace(/"/g, '&quot;')
            .replace(/'/g, '&#039;');
        };

        // Use the appropriate sanitizing function
        const sanitize = (window.htmlSanitizer && typeof window.htmlSanitizer.escapeHTML === 'function')
          ? window.htmlSanitizer.escapeHTML
          : escapeHtml;

        // If htmlSanitizer isn't available, log a warning once
        if (!window.htmlSanitizer && !window.warnedAboutMissingSanitizer) {
          console.warn('htmlSanitizer not available. Using basic string escaping for XSS protection.');
          window.warnedAboutMissingSanitizer = true;
        }

        // Return at most 5 results to avoid cluttering the command palette
        return sortedResults.slice(0, 5).map(result => {
          // Get raw values
          const rawTitle = typeof result.item.title === 'string' ? result.item.title : 'Untitled';
          const content = result.item.content || '';
          const rawExcerpt = result.item.excerpt || (content && content.substring(0, 100) + '...') || '';

          // Sanitize title and excerpt to prevent XSS
          const safeTitle = sanitize(rawTitle);
          const originalExcerpt = sanitize(rawExcerpt);

          return {
            id: `search-result-${result.refIndex}`,
            title: safeTitle,
            handler: () => {
              if (result.item.url) {
                window.location.href = result.item.url;
              }
            },
            section: "Search Results",
            icon: '<i class="fa-solid fa-file-lines"></i>',
            excerpt: originalExcerpt
          };
        });
      } catch (e) {
        console.error('Error performing search with Fuse:', e);
        return [];
      }
    }
  } catch (e) {
    console.error('Error searching database:', e);
    // Reset the fuseInitPromise if there was an error so future calls can try again
    window.searchHelper.fuseInitPromise = null;
  }

  return [];
};