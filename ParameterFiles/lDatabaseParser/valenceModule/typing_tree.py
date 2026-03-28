import os,sys


class tree_item:
#
#  tree item with these properties:
#  (1) level : set the ranking level in the tree
#  the larger the number is, the deeper the position is
#  (2) type_name : short_name of the type
#  (3) SMARTS : SMARTS string
#  (4) index : the index in the whole ranking tree
#  (5) local_index : the index in the same level connected to the same upper type
#  (6) upper_index : the upper level connected (global index)
#  (7) lower_index_dict : the lower level connected (global index)
#  (8) atom_type_index : the index of the atom type
#  (9) same_level_forward : the type forward in the same level
#  (10) same_level_backward : the type backward in the same level
#  (11) same_level_dict : the same level connected (global index)
#  (12) string : the string to describe the type
#
#  When adding the new item into the tree, the type_index and the type_local_index
#  will be automatically assigned based on your choice and the current status of the tree
#
#  Index follows the sequence of the lines in the type file, it is different from the atomic class number, because class
#  number could be shared by several lines which have different indices, i.e. their types are the same but the SMARTS and
#  description are different.
#

  def __init__(self, type_name, SMARTS, string, index, atom_type_index):
    self.level = 0
    self.name = type_name
    self.SMARTS = SMARTS
    self.index = index
    self.local_index = 0
    self.upper_index = -1
    self.lower_index_dict = {}
    self.lower_num = 0
    self.atom_type_index = atom_type_index
    self.same_level_forward = 0
    self.same_level_backward = 0
    self.same_level_dict = {}

class typing_tree:

  def __init__(self):
    self.level_num = 0
    self.type_num = 0
    self.premodified_items_set = {} # original items based on the natural sequence
    self.items_set = {} # items sorted by typing tree while the index is still based ontyping index from the typing file
    self.final_items_set = {} # items sorted by typing tree and the index has been changed based on the typing tree
    self.largest_level = 0

  def premodify_new_item(self, tree_item, type_level, upper_index):
    self.type_num += 1
    tree_item.level = type_level
    tree_item.upper_index = upper_index
    self.premodified_items_set[self.type_num] = tree_item

  def read_ranking_file(self, rank_file, type_file):
    f_0 = open(os.path.join(databasedir,rank_file))
    rank_lines = f_0.readlines()
    f_0.close()
    f_1 = open(os.path.join(databasedir,type_file))
    type_lines = f_1.readlines()
    f_1.close()
    items = {}

    for line in type_lines:
      terms = line.split()
      if(len(terms) == 0 or terms[0] == '#'):
        continue
      string = line.split('#')[-1]
      item = tree_item(terms[3], terms[0], string, int(terms[1]), int(terms[2]))
      items[int(terms[1])] = item

    for line in rank_lines:
      terms = line.split()
      type_level = int(terms[0])
      upper_index = int(terms[1])
      index = int(terms[2])
      self.premodify_new_item(items[index], type_level, upper_index)
      if(type_level > self.largest_level):
        self.largest_level = type_level

  def sorting_tree(self):
    for val in self.premodified_items_set.values():
      if(val.level not in self.items_set.keys()):
        self.items_set[val.level] = {}
      self.items_set[val.level][val.index] = val

    # sorted by typing tree and the index remains the initial style in typing file
    for i in range(2, self.largest_level+1):
      local_index = 0
      values = list(self.items_set[i].values())
      for j, item in enumerate(values):
        item.same_level_dict = self.items_set[i]
        if(j != 0):
          item.same_level_backward = values[j-1]
        if(j != len(values)-1):
          item.same_level_forward = values[j+1]
        local_index += 1
        item.local_index = local_index
        self.items_set[i-1][item.upper_index].lower_index_dict[item.index] = item

    # the index(keys of items_set) has been modified here
    total_index = 1
    for i in range(1, self.largest_level+1):
      dict_0 = {}
      for item in self.items_set[i].values():
        dict_0[total_index] = item
        total_index += 1
      self.final_items_set[i] = dict_0
    #print(self.items_set)
    #print(self.final_items_set)

  def same_level_search(self, searching_list, tree_item):
    forward_item = tree_item.same_level_forward
    backward_item = tree_item.same_level_backward
    local_index = tree_item.local_index
    if(forward_item == 0):
      for i in list(tree_item.same_level_dict.values())[:-1]:
        searching_list.append(i.atom_type_index)
    elif(backward_item == 0):
      for i in list(tree_item.same_level_dict.values())[1:]:
        searching_list.append(i.atom_type_index)
    else:
      for i in list(tree_item.same_level_dict.values())[local_index:]:
        searching_list.append(i.atom_type_index)
      for i in list(tree_item.same_level_dict.values())[:local_index-1]:
        searching_list.append(i.atom_type_index)

  def search_in_tree(self, atom_type_index):
    # return the searching list of the tree items
    find_item = 0 
    searching_list = []
    for item in self.premodified_items_set.items():
      if(item[0] == atom_type_index):
        find_item = item[1]
        break

    searching_list.append(find_item.atom_type_index)
    level = find_item.level
    current_item = find_item
    for i in range(level,1,-1):
      self.same_level_search(searching_list, current_item)
      current_item = self.premodified_items_set[current_item.upper_index]

    return searching_list


a = typing_tree()
global databasedir
moduledir = os.path.join(os.path.split(__file__)[0])
databasedir = moduledir.replace("valenceModule", 'dat')

a.read_ranking_file(os.path.join(databasedir, 'typing_tree.log'), os.path.join(databasedir, 'amoebaplusBondedType.dat'))
a.sorting_tree()
