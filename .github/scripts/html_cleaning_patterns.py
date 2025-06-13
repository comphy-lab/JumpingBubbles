#!/usr/bin/env python3
"""
HTML Cleaning Patterns Module

This module provides regex patterns for removing empty anchor tags from HTML files.
These empty anchors can cause JavaScript syntax errors in the generated documentation.

Usage:
    from html_cleaning_patterns import EMPTY_ANCHOR_PATTERNS, apply_empty_anchor_cleanup
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