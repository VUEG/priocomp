#!/usr/bin/python3
# -*- coding: utf-8 -*-


def pick_from_list(items, suffix):
    """ Pick an element from list ending with suffix.

    If no match is found, return None.

    :param items list of items to be searched.
    :suffix String suffix defining the match.
    """
    match = None
    for item in items:
        if item.endswith(suffix):
            match = item
    return match
