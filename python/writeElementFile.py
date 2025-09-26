# -*- coding: utf-8 -*-
"""
Created on Mon May  5 13:39:57 2025

@author: Tony
"""

def process_element_file(
    input_file: str,
    output_file: str,
    element_type: str,
    target_nodes: list,
    tag_tuple: tuple
):
    """
    Processes an element file to update the tag based on presence of certain nodes.

    Args:
        input_file (str): Path to the original .elem file.
        output_file (str): Path where the processed file will be saved.
        element_type (str): Type of element to filter on (e.g. "Tt").
        target_nodes (list): Node numbers to look for.
        tag_tuple (tuple): (default_tag, new_tag)
    """
    default_tag, new_tag = tag_tuple

    with open(input_file, 'r') as f:
        lines = f.readlines()

    # First line contains number of elements
    total_elements = int(lines[0].strip())
    updated_lines = [f"{total_elements}\n"]
    data_lines = lines[1:]

    # Convert target_nodes to a set for faster lookup
    target_nodes_set = set(target_nodes)

    for line in data_lines:
        parts = line.strip().split()
        if len(parts) != 6:
            continue  # skip malformed lines

        etype, *nodes_str, tag_str = parts
        if etype != element_type:
            updated_lines.append(line)
            continue

        nodes = list(map(int, nodes_str))
        tag = int(tag_str)

        # Update tag based on presence of target nodes
        if any(node in target_nodes_set for node in nodes):
            tag = new_tag
        else:
            tag = default_tag

        new_line = f"{etype} {' '.join(map(str, nodes))} {tag}\n"
        updated_lines.append(new_line)

    with open(output_file, 'w') as f_out:
        f_out.writelines(updated_lines)

    print(f"Processed file saved to: {output_file}")
