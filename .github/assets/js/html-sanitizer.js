/**
 * HTML Sanitizer Utility
 * Provides functions to safely escape HTML content to prevent XSS vulnerabilities
 */

// Create a namespace for our sanitizer to avoid global namespace pollution
window.htmlSanitizer = {
  /**
   * Escapes HTML special characters in a string to prevent XSS attacks
   * 
   * @param {string} str - The string to escape
   * @returns {string} - The escaped string
   */
  escapeHTML: function(str) {
    if (!str || typeof str !== 'string') {
      return '';
    }
    
    return str
      .replace(/&/g, '&amp;')
      .replace(/</g, '&lt;')
      .replace(/>/g, '&gt;')
      .replace(/"/g, '&quot;')
      .replace(/'/g, '&#039;');
  },
  
  /**
   * Creates a text node with the provided content to avoid HTML injection
   * 
   * @param {string} content - The content to create a text node for
   * @returns {Text} - A DOM Text node
   */
  createTextNode: function(content) {
    return document.createTextNode(content || '');
  },
  
  /**
   * Safely sets text content on an element
   * 
   * @param {HTMLElement} element - The element to set text on
   * @param {string} content - The text content to set
   */
  setTextContent: function(element, content) {
    if (!element || !element.nodeType) {
      return;
    }
    
    // Clear any existing content
    while (element.firstChild) {
      element.removeChild(element.firstChild);
    }
    
    // Add the new content as a text node
    element.appendChild(this.createTextNode(content));
  }
};