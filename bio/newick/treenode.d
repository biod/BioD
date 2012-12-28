module bio.newick.treenode;

/** Node in a Newick tree */
struct NewickNode {

    /** Parent node (null if node is root) */
    NewickNode* parent;

    /** Node name (null if node is not named) */
    string name;

    /** Distance to parent node */
    double distance_to_parent = 1.0;

    /** Child nodes */
    NewickNode*[] children;
}
