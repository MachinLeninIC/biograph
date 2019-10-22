
import tempfile
import subprocess
import os
import Bio.PDB as PDB

class CDHitGroup:

    @staticmethod
    def get_group(proteins, similarity=0.9, word_size=5, memory_mb=1024, threads=1):
        """Given a list of pyprot.Protein objects, runs CDHit
        to cluster sequences based on similarity.
        CDHit must be installed for this method to work."""
        builder = PDB.Polypeptide.PPBuilder()
        sequences = []
        names = []
        for p in proteins:
            chains = p.pdb[0].get_list()
            for i, pp in enumerate(builder.build_peptides(p.pdb, aa_only=0)):
                sequences.append(str(pp.get_sequence()))
                try:
                    names.append("{}_{}".format(p.pdb.id, chains[i].id))
                except IndexError:
                    print("Index error while processing protein ID {}".format(p.pdb.id))
                    print("Check that biopython finds as many sequences as chains there are")

        temppath = tempfile.mkdtemp()
        seqfile = os.path.join(temppath, "sequences.fasta.txt")
        with open(seqfile, "w") as f:
            for name, seq in zip(names, sequences):
                f.write(">{}\n{}\n".format(name, seq))

        #TODO: check where cdhit is, or assert that cdhit is installed.
        out = subprocess.run(["./cdhit/cd-hit",
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

        groups = [grouping[name] for name in names]
        return groups