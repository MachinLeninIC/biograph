
import tempfile
import subprocess
import os
import Bio.PDB as PDB
from biograph.UnionFind import UnionFind

class CDHitGroup:

    @staticmethod
    def get_protein_groups(protein_chain_map, chain_group_map):
        """Given chain-level groups and a protein-chain map, create
        groups so that no two proteins that share a chain group are
        in different protein groups.
        Parameters
        ----------
        protein_chain_map: dict str => list
            Dictionary mapping proteins to a list of chain names.
        chain_group_map: dict str => int
            Dictionary mapping chain names to group ids.
        """
        # UnionFind uses numbers 0..n-1 to identify groups, so map
        # cluster ids to that.
        group_id = {}
        inverse_group_id = {}
        for group in chain_group_map.values():
            group_id[group] = group_id.get(group, len(group_id))
            inverse_group_id[group_id[group]] = group

        uf = UnionFind(len(group_id))
        for protein, chains in protein_chain_map.items():
            first_group = chain_group_map[chains[0]]
            for chain in chains:
                chain_group = chain_group_map[chain]
                uf.union(group_id[first_group], group_id[chain_group])

        # Now that UF group ids are stable we can assign protein groups.
        protein_group_map = {}
        for protein, chains in protein_chain_map.items():
            uf_group = uf.find(group_id[chain_group_map[chains[0]]])
            protein_group_map[protein] = inverse_group_id[uf_group]

        return protein_group_map

    @staticmethod
    def get_group_by_sequences(sequences, names, similarity=0.9, word_size=5, memory_mb=1024, threads=1):
        """Same as get_group but manually providing sequences and names."""
        temppath = tempfile.mkdtemp()
        seqfile = os.path.join(temppath, "sequences.fasta.txt")
        with open(seqfile, "w") as f:
            for name, seq in zip(names, sequences):
                f.write(">{}\n{}\n".format(name, seq))

        MODULEDIR = os.path.dirname(os.path.abspath(__file__))
        print("Using MODULEDIR="+MODULEDIR)
        out = subprocess.run([
            os.path.join(MODULEDIR, "cdhit", "cd-hit"),
            "-i", seqfile, "-o", "output.txt",
            "-c", str(similarity), "-n", str(word_size),
            "-M", str(memory_mb), "-T", str(threads),
            "-bak", "1", "-d", "0"])

        if out.returncode != 0:
            print("Error running CDHit")
            return

        grouping = {}
        with open("output.txt.clstr") as f:
            """Lines are of the format:
            >Cluster 2
            0       164aa, >6std_ASDFGHJ... *
            1       164aa, >6std_B... at 100.00%"""
            line = f.readline()
            current_cluster = -1
            while line:
                if line[0] == ">":
                    current_cluster += 1
                else:
                    name = line[line.find(">")+1:line.find("...")]
                    grouping[name] = current_cluster
                line = f.readline()

        # We use .get() since some chains passed to CDHit may not be in
        # any cluster. This can happen if a chain is just a peptide or
        # if it is empty. .get() returns None if it is not found, so
        # the resulting list has the same len() as names.
        groups = [grouping.get(name) for name in names]
        return groups
    @staticmethod
    def get_group(proteins, similarity=0.9, word_size=5, memory_mb=1024, threads=1):
        """Given a list of biograph.Protein objects, runs CDHit
        to cluster sequences based on similarity.
        CDHit must be installed for this method to work."""
        sequences = []
        names = []
        for p in proteins:
            sequences.extend(p.sequences.values())
            names.extend(p.sequences.keys())

        return CDHitGroup.get_group_by_sequences(sequences, names, similarity=similarity,
            word_size=word_size, memory_mb=memory_mb, threads=threads)