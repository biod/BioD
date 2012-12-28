module bio.newick.treenode;

/** Node in a Newick tree */
struct NewickNode {

    /** Node name (null if node is not named) */
    string name;

    /** Distance to parent node (NaN if not specified) */
    double distance_to_parent;

    /** Child nodes */
    NewickNode*[] children;
}
