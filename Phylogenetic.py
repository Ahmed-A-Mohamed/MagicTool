# you have to use the following commends to use these package
# you have to install matpotlib
# you have to install biopythone
# sudo apt-get install python3-tk



dnd = input("Enter .dnd")


from Bio import Phylo

tree = Phylo.read(dnd, "newick")
Phylo.draw_ascii(tree)
tree.rooted = True
tree.root.color = "salmon"
tree.clade[0, 1].color = "blue"
Phylo.draw(tree, branch_labels=lambda c: c.branch_length)

