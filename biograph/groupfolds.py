
import tempfile
import subprocess
import os
import Bio.PDB as PDB

class CDHitGroup:

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