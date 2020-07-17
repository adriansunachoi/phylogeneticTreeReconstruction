"""
Python Tree Class
Stuart Bradley - 5931269
23-05-2014
"""

class Tree:

	def __init__(self, root=None):
		self.root = root 

	def set_root(self, root):
		self.root = root

	def get_root(self):
		return self.root

	def get_leaves(self):
		return self.root.get_leaves()

	def get_leaf_count(self):
		return len(self.get_leaves())
    
   	def get_newick(self):
   		return self.root.get_newick() + ";"

   	def __str__(self):
   		return self.get_newick()