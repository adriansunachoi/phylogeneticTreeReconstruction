""" 
Python Node Class
Stuart Bradley - 5931269
23-05-2014
"""

class Node:
	def __init__(self, label=None):
		self.parent = None
		self.children = []
		self.height = -1.0
		self.label = label
		self.sequence = None

	def get_parent(self):
		return self.parent

	def set_parent(self, parent):
		self.parent = parent

	def get_children(self):
		return self.children

	def add_child(self, child):
		self.children.append(child)
		child.set_parent(self)

	def remove_child(self, child):
		self.children.remove(child)

	def set_height(self, height):
		self.height = height

	def get_height(self):
		return self.height

	def is_root(self):
		return self.parent == None

	def is_leaf(self):
		return not self.children

	def get_sequence(self):
		return self.sequence

	def set_sequence(self, sequence):
		self.sequence = sequence

	def get_label(self):
		return self.label

	def set_label(self, label):
		self.label = label

	def get_leaves(self):
		leaf_list = []

		if (self.is_leaf()):
			leaf_list.append(self)
		else:
			for child in self.children:
				leaf_list.extend(child.get_leaves())

		return leaf_list

	def get_newick(self):
		sb = ""

		if (not self.is_leaf()):
			sb += "("
			for i in range(0, len(self.children)):
				if (i>0):
					sb += ","
				sb += self.children[i].get_newick()
			sb += ")"
		
		if (self.label != None):
			sb += self.label

		branch_length = -1.0
		if (not self.is_root()):
			branch_length = self.parent.height - self.height
		else:
			branch_length = 0.0

		sb += ":" + str(branch_length)

		return sb