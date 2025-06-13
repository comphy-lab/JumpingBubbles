/* ===================================================================
 * Main JS
 * ------------------------------------------------------------------- */

(function(html) {

    'use strict';
    
    // Use centralized DEBUG flag from config
    const DEBUG = window.CONFIG ? window.CONFIG.DEBUG : false;
    
    /**
     * Copies the specified text to the clipboard, using the Clipboard API if available, with a fallback for older browsers.
     *
     * Updates the button state to indicate success or failure.
     *
     * @param {string} text - The text to be copied to the clipboard.
     * @param {HTMLElement} button - The button element to update based on the copy result.
     */
    function copyToClipboard(text, button) {
        // Try to use modern Clipboard API first with optional chaining
        if (navigator.clipboard?.writeText) {
            navigator.clipboard.writeText(text)
                .then(() => updateButtonState(button, true))
                .catch(err => {
                    if (DEBUG) {
                        console.error('Clipboard API failed:', err);
                    }
                    fallbackCopyToClipboard(text, button);
                });
        } else {
            // Fall back to deprecated execCommand method for older browsers
            fallbackCopyToClipboard(text, button);
        }
    }
    
    /**
     * Attempts to copy text to the clipboard using a hidden textarea and the deprecated `execCommand('copy')` method.
     *
     * Updates the button state to indicate success or failure.
     *
     * @param {string} text - The text to be copied to the clipboard.
     * @param {HTMLElement} button - The button element to update based on the copy result.
     *
     * @remark This function is used as a fallback when the modern Clipboard API is unavailable.
     */
    function fallbackCopyToClipboard(text, button) {
        const textarea = document.createElement('textarea');
        textarea.value = text;
        textarea.style.position = 'fixed';
        textarea.style.left = '-9999px';
        // Add to body temporarily
        document.body.appendChild(textarea);
        
        try {
            textarea.select();
            const success = document.execCommand('copy');
            if (success) {
                updateButtonState(button, true);
            } else {
                if (DEBUG) {
                    console.error('Copy command failed');
                }
                updateButtonState(button, false);
            }
        } catch (err) {
            if (DEBUG) {
                console.error('Fallback copy failed:', err);
            }
            updateButtonState(button, false);
        } finally {
            document.body.removeChild(textarea);
        }
    }
    
    function updateButtonState(button, success) {
        if (success) {
            const icon = button.querySelector('i');
            button.classList.add('copied');
            
            if (icon) {
                icon.classList.remove('fa-copy');
                icon.classList.add('fa-check');
            }
            
            // Reset button state after 2 seconds
            setTimeout(() => {
                button.classList.remove('copied');
                if (icon) {
                    icon.classList.remove('fa-check');
                    icon.classList.add('fa-copy');
                }
            }, 2000);
        }
    }

    /* Preloader
    * -------------------------------------------------- */
    const preloader = document.querySelector("#preloader");
    if (preloader) {
        window.addEventListener('load', function() {
            document.querySelector('body').classList.remove('ss-preload');
            document.querySelector('body').classList.add('ss-loaded');
            preloader.style.display = 'none';
        });
    }
    
    // No need for a resize event handler as the CSS will handle everything

    // Only load content if the functions exist
    if (typeof loadAboutContent === 'function') {
        window.addEventListener('load', loadAboutContent);
    }
    if (typeof loadNewsContent === 'function') {
        window.addEventListener('load', loadNewsContent);
    }

    /* Mobile Menu
    * -------------------------------------------------- */
    const menuToggle = document.querySelector('.s-header__menu-toggle');
    const nav = document.querySelector('.s-header__nav');
    const closeBtn = document.querySelector('.s-header__nav-close-btn');
    const menuLinks = document.querySelectorAll('.s-header__nav-list a');
    
    // Debug elements
    if (DEBUG) {
        console.log('Menu elements found:', { 
            menuToggle: menuToggle !== null, 
            nav: nav !== null, 
            closeBtn: closeBtn !== null,
            menuLinksCount: menuLinks ? menuLinks.length : 0
        });
    }

    // Handle click outside
    document.addEventListener('click', function(e) {
        if (nav && nav.classList.contains('is-active')) {
            // Check if click is outside nav and not on menu toggle
            if (!nav.contains(e.target) && !menuToggle?.contains(e.target)) {
                if (DEBUG) {
                    console.log('Click outside detected');
                }
                nav.classList.remove('is-active');
                // Reset the style
                nav.style.right = '-300px';
            }
        }
    });

    if (menuToggle) {
        menuToggle.addEventListener('click', function(e) {
            e.preventDefault();
            e.stopPropagation(); // Prevent document click from immediately closing
            if (DEBUG) {
                console.log('Menu toggle clicked');
            }
            if (nav) {
                if (DEBUG) {
                    console.log('Adding is-active class to nav');
                }
                nav.classList.add('is-active');
                
                // Make sure the style change is applied
                nav.style.right = '0';
            } else {
                if (DEBUG) {
                    console.error('Nav element not found when toggle clicked');
                }
            }
        });
    } else {
        if (DEBUG) {
            console.error('Menu toggle element not found');
        }
    }

    if (closeBtn) {
        closeBtn.addEventListener('click', function(e) {
            e.preventDefault();
            if (DEBUG) {
                console.log('Close button clicked');
            }
            if (nav) {
                if (DEBUG) {
                    console.log('Removing is-active class from nav');
                }
                nav.classList.remove('is-active');
                // Reset the style
                nav.style.right = '-300px';
            } else {
                if (DEBUG) {
                    console.error('Nav element not found when close clicked');
                }
            }
        });
    } else {
        if (DEBUG) {
            console.error('Close button element not found');
        }
    }

    if (menuLinks && menuLinks.length > 0) {
        menuLinks.forEach(link => {
            link.addEventListener('click', () => {
                if (DEBUG) {
                    console.log('Menu link clicked');
                }
                if (nav) {
                    if (DEBUG) {
                        console.log('Removing is-active class from nav');
                    }
                    nav.classList.remove('is-active');
                    // Reset the style
                    nav.style.right = '-300px';
                }
            });
        });
    }

    /* Smooth Scrolling
    * -------------------------------------------------- */
    document.querySelectorAll('a[href^="#"]').forEach(anchor => {
        anchor.addEventListener('click', function (e) {
            e.preventDefault();
            const href = this.getAttribute('href');
            
            // Skip navigation placeholders like #0
            if (href === '#0' || href === '#') {
                return;
            }
            
            const target = document.querySelector(href);
            if (target) {
                target.scrollIntoView({
                    behavior: 'smooth'
                });
            }
        });
    });

    /* Back to Top
    * -------------------------------------------------- */
    const goTop = document.querySelector('.ss-go-top');

    if (goTop) {
        window.addEventListener('scroll', function() {
            if (window.pageYOffset > 800) {
                goTop.classList.add('link-is-visible');
            } else {
                goTop.classList.remove('link-is-visible');
            }
        });
    }

    document.addEventListener('DOMContentLoaded', function() {
        const images = document.querySelectorAll('.member-image img[loading="lazy"]');
        
        images.forEach(img => {
            if (img.complete) {
                img.parentElement.classList.add('loaded');
            } else {
                img.addEventListener('load', function() {
                    img.parentElement.classList.add('loaded');
                });
            }
        });

        // Email copy functionality
        const copyButtons = document.querySelectorAll('.copy-btn');
        copyButtons.forEach(button => {
            button.addEventListener('click', function() {
                const textToCopy = this.getAttribute('data-clipboard-text');
                copyToClipboard(textToCopy, this);
            });
        });

        // Add accessible names to all copy buttons on document load
        copyButtons.forEach(button => {
            // Get the email text from data-text or data-clipboard-text attribute
            const emailText = button.getAttribute('data-text') || button.getAttribute('data-clipboard-text');
            // Add aria-label if it doesn't exist
            if (!button.hasAttribute('aria-label') && emailText) {
                button.setAttribute('aria-label', `Copy email address ${emailText}`);
            }
        });
    });

    /* Copy Email Functionality
    * -------------------------------------------------- */
    window.copyEmail = function(button) {
        const text = button.getAttribute('data-text') || button.getAttribute('data-clipboard-text');
        copyToClipboard(text, button);
    };

})(document.documentElement);