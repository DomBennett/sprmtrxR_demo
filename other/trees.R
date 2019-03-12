# Download the primates tree and reformat
# Perelman et al. 2011
# https://treebase.org/treebase-web/search/study/trees.html?id=12186
primates_tree <- ape::read.nexus(file.path('0_data', 'T50066.nex'))
match <- regexpr(pattern = '^[a-zA-Z]+_[a-zA-Z]+',
                 text = primates_tree$tip.label)
match_lengths <- attr(match, 'match.length')
new_tip_labels <- substr(x = primates_tree$tip.label, start = 1,
                         stop = match_lengths)
to_drop <- primates_tree$tip.label[duplicated(new_tip_labels)]
primates_tree <- ape::drop.tip(primates_tree, tip = to_drop)
new_tip_labels <- new_tip_labels[!duplicated(new_tip_labels)]
primates_tree$tip.label <- new_tip_labels
ape::write.tree(phy = primates_tree, file = file.path('0_data', 'primate.tre'))
