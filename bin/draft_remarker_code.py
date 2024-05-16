import re


def extract_mask_coordinates(mask):
    start, length = map(int, mask.strip('N').split('+'))
    return list(range(start, start + length))


def get_mask_coordinates(pseudo_cigar):
    # Get masking coordinates
    masking_pattern = r'\d+\+\d+N'
    masks = re.findall(masking_pattern, pseudo_cigar)
    mask_coordinates = []
    for mask in masks:
        mask_coordinates.extend(extract_mask_coordinates(mask))
    return set(mask_coordinates)


def extract_codon_from_asv(reference, sequence, mask_coordinates, start_position):
    preceding_indel_present = False

    # No indels
    if '-' not in reference and '-' not in sequence:
        codon = sequence[start_position:start_position+3]
    # Indels
    else:
        index = 0
        while index < start_position:
            ref_base = reference[index]
            seq_base = sequence[index]
            if ref_base == '-':
                start_position += 1
                if index + 1 not in mask_coordinates:
                    preceding_indel_present = True
            if seq_base == '-':
                if index + 1 not in mask_coordinates:
                    preceding_indel_present = True
            index += 1
        codon = sequence[start_position:start_position+3]
    return codon, preceding_indel_present


def process_resmarker(reference, sequence, start_position, pseudo_cigar, ref_codon):
    if pseudo_cigar == '.':
        codon_masked = False
        codon = ref_codon
        preceding_indel_present = False
    else:
        codon_coordinates = list(range(start_position, start_position+3))
        # Need to change this
        mask_coordinates = get_mask_coordinates(pseudo_cigar)

        # Check if there's any overlap
        codon_masked = bool(mask_coordinates.intersection(codon_coordinates))
        # Get codon from asv
        codon, preceding_indel_present = extract_codon_from_asv(
            reference, sequence, mask_coordinates, start_position)
    return codon, codon_masked, preceding_indel_present


# REF
reference = 'ACTGACTAAG'
sequence = 'ACTGACTAAG'
start_position = 6
print('REF')
print(process_resmarker(reference, sequence,
      start_position, '.', 'TAA'))  # Output: TGA
print('')

# NO INDEL
reference = 'ACTGACTAAG'
sequence = 'ACTGACTGAG'
start_position = 6
print('No indel')
print(process_resmarker(reference, sequence,
      start_position, '8G', 'TAA'))  # Output: TGA
print('')

# INSERTION outside of mask - move coordinates by + one
reference = 'ACTGA-CTAAG'
sequence = 'ACTGAACTGAG'
start_position = 6
print('Insertion outside of mask')
print(process_resmarker(reference, sequence,
      start_position, '6I=A8G', 'TAA'))  # Output: TGA
print('')

# DELETION outside of mask - don't move coordinates
reference = 'ACTGACTAAG'
sequence = 'ACTG-CTGAG'
start_position = 6
print('Deletion outside of mask')
print(process_resmarker(reference, sequence,
      start_position, '5D=A8G', 'TAA'))  # Output: TGA
print('')

# INSERTION after codon
reference = 'ACTGACTAA-G'
sequence = 'ACTGACTGATG'
start_position = 6
print('Insertion after codon')
print(process_resmarker(reference, sequence,
      start_position, '6I=A8G', 'TAA'))  # Output: TGA
print('')

# DELETION after codon
reference = 'ACTGACTAAG'
sequence = 'ACTGACTGA-'
start_position = 6
print('Deletion after codon')
print(process_resmarker(reference, sequence,
      start_position, '5D=A8G', 'TAA'))  # Output: TGA
print('')

# MASKED INSERTION AND PARTIAL CODON
reference = 'ACTGA-CTAAG'
sequence = 'ACTGAACTGAG'
start_position = 6
print('Masked insertion and masked partial codon')
print(process_resmarker(reference, sequence,
      start_position, '5+3N8G', 'TAA'))  # Output: TGA
print('')

# MASKED DELETION AND CODON
reference = 'ACTGACTAAG'
sequence = 'ACTG-CTGAG'
start_position = 6
print('Masked deletion and masked codon')
print(process_resmarker(reference, sequence,
      start_position, '4+3N8G', 'TAA'))  # Output: TGA
print('')

# Insertion in the codon
reference = 'ACTGACT-AAG'
sequence = 'ACTGACTAGAG'
start_position = 6
print('Insertion in the codon')
print(process_resmarker(reference, sequence,
      start_position, '5+3N8I=A8G', 'TAA'))  # Output: TAG
print('')

# Deletion in the codon
reference = 'ACTGACTAAG'
sequence = 'ACTGACT-AG'
start_position = 6
print('Deletion in the codon')
print(process_resmarker(reference, sequence,
      start_position, '5+3N8D=A', 'TAA'))  # Output: T-A
print('')

# COMBO - masked deletion, unmasked insertion, masked codon
reference = 'ACT--GACTAAG'
sequence = '-CTTTGACT-AG'
start_position = 6
print('Deletion in the codon')
print(process_resmarker(reference, sequence,
      start_position, '1+2N4I=28G7+4N', 'TAA'))  # Output: T-A

# Need to handle asv not being long enough to get codon
