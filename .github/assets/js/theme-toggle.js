// Theme toggle functionality
document.addEventListener('DOMContentLoaded', () => {
  // Check for saved theme preference or use OS preference
  const getThemePreference = () => {
    // Check if user has already selected a theme
    try {
      const savedTheme = localStorage.getItem('theme');
      if (savedTheme) {
        return savedTheme;
      }
    } catch (e) {
      // localStorage access failed, fall back to OS preference
      console.warn('localStorage access failed:', e);
    }
    // If not saved or localStorage failed, check OS preference
    return window.matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light';
  };
  
  // Set theme on document
  const setTheme = (theme) => {
    document.documentElement.setAttribute('data-theme', theme);
    try {
      localStorage.setItem('theme', theme);
    } catch (e) {
      // localStorage access failed, theme will not persist
      console.warn('Failed to save theme preference:', e);
    }
  };
  
  // Initialize theme
  setTheme(getThemePreference());
  
  // Add click event listener to theme toggle
  const themeToggle = document.getElementById('theme-toggle');
  if (themeToggle) {
    themeToggle.addEventListener('click', () => {
      const currentTheme = document.documentElement.getAttribute('data-theme') || 'light';
      const newTheme = currentTheme === 'light' ? 'dark' : 'light';
      setTheme(newTheme);
    });
  }
});
