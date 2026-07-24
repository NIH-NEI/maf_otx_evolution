import requests
import sys

def read_dna_sequence(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()
    # Skip FASTA header if present
    dna = ''.join(line.strip().replace(' ', '').replace('-', '') for line in lines if not line.startswith('>'))
    print(dna)
    return dna.upper()

def reverse_complement(seq):
    comp = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
    return ''.join(comp.get(base, 'N') for base in reversed(seq))

def search_ucsc_blat(dna_seq, genome='mm9'):
    url = 'https://genome.ucsc.edu/cgi-bin/hgBlat'
    params = {
        'userSeq': dna_seq,
        'db': genome,
        'type': 'DNA',
        'output': 'json',
        'Submit': 'submit'
    }

    response = requests.post(url, data=params)
    if response.status_code == 200:
        if "no matches" in response.text.lower():
            return None
        else:
            return response.text
    else:
        print(f"Error contacting UCSC BLAT server: {response.status_code}")
        return None

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python search_dna_ucsc.py <query_file.txt> <ucsc_genome>")
        sys.exit(1)

    query_file = sys.argv[1]
    genome = sys.argv[2]
    seq = read_dna_sequence(query_file)
    print("Query length:", len(seq))

    print("Searching forward strand...")
    result_fwd = search_ucsc_blat(seq, genome = genome)
    if result_fwd:
        print("\n🔍 Match found on forward strand:")
        print(result_fwd)
    else:
        print("No match found on forward strand.")

    print("\nSearching reverse strand...")
    rev_seq = reverse_complement(seq)
    result_rev = search_ucsc_blat(rev_seq)
    if result_rev:
        print("\n🔍 Match found on reverse strand:")
        print(result_rev)
    else:
        print("No match found on reverse strand.")
