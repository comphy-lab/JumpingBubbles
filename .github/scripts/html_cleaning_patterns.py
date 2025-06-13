#!/usr/bin/env python3
"""
HTML Cleaning Patterns Module

This module provides a collection of regex patterns for HTML cleaning tasks,
ensuring consistent behavior across different HTML cleaning scripts.

Usage:
    from html_cleaning_patterns import EMPTY_ANCHOR_PATTERNS, SANITIZE_PATTERNS
"""

# Patterns for identifying and removing empty anchor tags
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

# Patterns for sanitizing HTML content to prevent XSS
SANITIZE_PATTERNS = {
    # Script tag patterns
    'script_open': r'<\s*script',
    'script_close': r'<\s*\/\s*script',
    
    # Iframe tag patterns
    'iframe_open': r'<\s*iframe',
    'iframe_close': r'<\s*\/\s*iframe',
    
    # Event handler attributes
    'event_handlers': r'on\w+\s*=\s*["\'][^"\']*["\']',
    
    # JavaScript URLs
    'javascript_urls': r'javascript\s*:'
}

# Replacement values for sanitization
SANITIZE_REPLACEMENTS = {
    'script_open': '&lt;script',
    'script_close': '&lt;/script',
    'iframe_open': '&lt;iframe',
    'iframe_close': '&lt;/iframe',
    'event_handlers': '',
    'javascript_urls': 'disabled-javascript:'
}

# Utility functions for common operations
def apply_empty_anchor_cleanup(content):
    """
    Apply all empty anchor cleanup patterns to the content.
    
    Args:
        content: HTML content to clean
        
    Returns:
        str: Cleaned HTML content with empty anchors removed
    """
    import re
    result = content
    for pattern in EMPTY_ANCHOR_PATTERNS:
        result = re.sub(pattern, '', result)
    return result

def apply_html_sanitization(content):
    """
    Apply all HTML sanitization patterns to the content.
    
    Args:
        content: HTML content to sanitize
        
    Returns:
        str: Sanitized HTML content
    """
    import re
    import html
    
    result = content
    for key, pattern in SANITIZE_PATTERNS.items():
        replacement = SANITIZE_REPLACEMENTS[key]
        result = re.sub(pattern, replacement, result, flags=re.IGNORECASE)
    
    # Unescape first to prevent double-encoding from previous regex replacements
    result = html.unescape(result)
    
    # Finally, escape any remaining HTML special characters
    result = html.escape(result, quote=False)  # Don't escape quotes again
    
    return result