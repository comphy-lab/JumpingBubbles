#!/usr/bin/env python3
"""
HTML Cleaning Patterns Module

This module provides regex patterns for removing empty anchor tags from HTML files.
These empty anchors can cause JavaScript syntax errors in the generated documentation.

Usage:
    from html_cleaning_patterns import EMPTY_ANCHOR_PATTERNS, apply_empty_anchor_cleanup
"""

import re

# Pre-compiled patterns for identifying and removing empty anchor tags
# These are compiled at import time for better performance when processing multiple files
EMPTY_ANCHOR_PATTERNS = [
    # Combined pattern for id and href attributes in either order with optional newlines
    re.compile(r'<a\s+id=[\'"]?([^\s>]*)[\'"]?\s+href=[\'"]?#[\'"]?\s*>\s*(?:\n\s*)*</a>', re.IGNORECASE),
    re.compile(r'<a\s+href=[\'"]?#[\'"]?\s+id=[\'"]?([^\s>]*)[\'"]?\s*>\s*(?:\n\s*)*</a>', re.IGNORECASE),
    
    # Single pattern for id attribute only with optional newlines
    re.compile(r'<a\s+id=[\'"]?([^\s>]*)[\'"]?\s*>\s*(?:\n\s*)*</a>', re.IGNORECASE),
    
    # Single pattern for href='#' only with optional newlines
    re.compile(r'<a\s+href=[\'"]?#[\'"]?\s*>\s*(?:\n\s*)*</a>', re.IGNORECASE),
    
    # Combined pattern for unquoted attributes in either order
    re.compile(r'<a\s+(?:id=([^\s>]*)\s+href=#|href=#\s+id=([^\s>]*))\s*>\s*(?:\n\s*)*</a>', re.IGNORECASE)
]

def apply_empty_anchor_cleanup(content):
    """
    Apply all empty anchor cleanup patterns to the content.
    
    Args:
        content: HTML content to clean
        
    Returns:
        str: Cleaned HTML content with empty anchors removed
    """
    result = content
    for pattern in EMPTY_ANCHOR_PATTERNS:
        result = pattern.sub('', result)
    return result