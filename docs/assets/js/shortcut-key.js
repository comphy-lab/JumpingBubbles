document.addEventListener('DOMContentLoaded', function() {
  // Detect if the user is on a Mac using modern API with fallback
  let isMac = false;
  
  if (navigator.userAgentData && navigator.userAgentData.platform) {
    // Use modern User-Agent Client Hints API when available
    isMac = navigator.userAgentData.platform.toUpperCase().indexOf('MAC') >= 0;
  } else {
    // Fall back to the deprecated navigator.platform for older browsers
    isMac = navigator.platform.toUpperCase().indexOf('MAC') >= 0;
  }
  
  // Update the displayed shortcut text
  const defaultThemeElements = document.querySelectorAll('.default-theme-text');
  const macThemeElements = document.querySelectorAll('.mac-theme-text');
  
  defaultThemeElements.forEach(function(element) {
    element.style.display = isMac ? 'none' : 'inline';
  });
  
  macThemeElements.forEach(function(element) {
    element.style.display = isMac ? 'inline' : 'none';
  });
}); 