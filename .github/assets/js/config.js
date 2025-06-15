/**
 * Global configuration settings for the website
 * 
 * This file centralizes configuration settings that are used across
 * multiple JavaScript files to improve maintainability and consistency.
 */

(function() {
  // Get the current environment - defaults to production
  const isDevEnvironment = function() {
    // Check if we're in a development environment
    // This can be expanded with more checks as needed
    const hostname = window.location.hostname;
    return hostname === 'localhost' || 
           hostname === '127.0.0.1' || 
           hostname.includes('192.168.') ||
           hostname.includes('.local') ||
           hostname.includes('gitpod.io') ||
           hostname.includes('github.dev');
  };

  // Global configuration object
  window.CONFIG = {
    // Set DEBUG based on environment
    DEBUG: isDevEnvironment(),
    
    // Other global configuration options can be added here
    VERSION: '1.0.0',
    
    // Utility function to access debug state with logging disabled for production
    isDebug: function() {
      return this.DEBUG;
    },
    
    // Debug logging function that respects DEBUG setting
    log: function(message, data) {
      if (this.DEBUG) {
        if (data !== undefined) {
          console.log(message, data);
        } else {
          console.log(message);
        }
      }
    },
    
    // Error logging always happens regardless of DEBUG
    error: function(message, error) {
      if (error !== undefined) {
        console.error(message, error);
      } else {
        console.error(message);
      }
    }
  };
})();