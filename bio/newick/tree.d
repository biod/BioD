module bio.newick.tree;

import bio.newick.parser;
import bio.newick.treenode;

/** 
  Newick tree
*/
class NewickTree {

    private {
        NewickNode* _root;
    }

    /** Root node, may be null */
    NewickNode* root() @property {
        return _root;
    }

    /** Creates a tree from its string representation */
    static NewickTree fromString(string str) {
        auto tree = new NewickTree();
        tree._root = bio.newick.parser.parse(str); 
        tree._fillDictionary();
        return tree;
    }

    /** Creates a tree with a given root node */
    static NewickTree fromRootNode(NewickNode* root) {
        auto tree = new NewickTree();
        tree._root = root;
        tree._fillDictionary();
        return tree;
    }

    /** Dictionary of nodes */
    NewickNode*[string] nodes;

    /** Distance between two nodes */
    double distance(NewickNode* n1, NewickNode* n2) {
        if (n1 is null || n2 is null) return double.nan;
        if (n1 == n2) return 0;
        if (n1 == root) return distance(root, n2.parent) + n2.distance_to_parent;
        if (n2 == root) return distance(n1.parent, root) + n1.distance_to_parent;
        assert(n1.parent !is null);
        assert(n2.parent !is null);
        return distance(n1.parent, n2.parent) + n1.distance_to_parent 
                                              + n2.distance_to_parent;
    }

    /** ditto */
    double distance(string n1, string n2) {
        return distance(nodes[n1], nodes[n2]);
    }

    private {
        void _fillDictionary() {
            _fillDictRecur(root);
        }

        void _fillDictRecur(NewickNode* node) {
            if (node !is null) {
                if (node.name !is null)
                    nodes[node.name] = node;

                foreach (child; node.children) {
                    _fillDictRecur(child);
                }
            }
        }
    }
}

unittest {
    import std.math;

    auto tree = NewickTree.fromString("(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);");
    assert(tree.root !is null);
    assert(tree.root.name is null);
    assert(approxEqual(tree.distance("A", "B"), 0.3));
    assert(approxEqual(tree.distance("B", "C"), 1.0));
    assert(approxEqual(tree.distance("C", "D"), 0.7));
    assert(approxEqual(tree.distance("A", "C"), 0.9));
    assert(approxEqual(tree.distance("B", "D"), 1.1));
    assert(approxEqual(tree.distance("A", "D"), 1.0));
    assert(tree.root.children.length == 3);
    assert(tree.root.children[0].name == "A");
    assert(approxEqual(tree.root.children[0].distance_to_parent, 0.1));
    assert(tree.root.children[1].name == "B");
    assert(approxEqual(tree.root.children[1].distance_to_parent, 0.2));
    auto child = tree.root.children[2];
    assert(child !is null);
    assert(child.parent == tree.root);
    assert(child.name is null);
    assert(approxEqual(child.distance_to_parent, 0.5));
    assert(child.children.length == 2);
    assert(child.children[0].name == "C");
    assert(approxEqual(child.children[0].distance_to_parent, 0.3));
    assert(child.children[1].name == "D");
    assert(approxEqual(child.children[1].distance_to_parent, 0.4));

    tree = NewickTree.fromString("(A,B,(C,D));");
    assert(tree.root !is null);
    assert(tree.root.children.length == 3);
    assert(tree.root.children[0].name == "A");
    assert(tree.root.children[1].name == "B");
    assert(tree.root.children[2].name is null);
    assert(tree.root.children[2].children[0].name == "C");
    assert(tree.root.children[2].children[1].name == "D");
    assert(tree.nodes["D"].parent.children[0].name == "C");

    tree = NewickTree.fromString("((B:0.2,(C:0.3,D:0.4)E:0.5)F:0.1)A;");
    assert(tree.root.name == "A"); 
    assert(tree.root.children.length == 1);
    child = tree.root.children[0];
    assert(child.parent == tree.root);
    assert(child.name == "F");
    assert(approxEqual(child.distance_to_parent, 0.1));
    assert(child.children.length == 2);
    assert(child.children[0].name == "B");
    assert(approxEqual(child.children[0].distance_to_parent, 0.2));
    child = child.children[1];
    assert(child.name == "E");
    assert(approxEqual(child.distance_to_parent, 0.5));
    assert(child.children.length == 2);
    assert(child.children[0].name == "C");
    assert(approxEqual(child.children[0].distance_to_parent, 0.3));
    assert(child.children[1].name == "D");
    assert(approxEqual(child.children[1].distance_to_parent, 0.4));

    tree = NewickTree.fromString("(,,(,));");
    assert(tree.root !is null);
    assert(tree.root.children.length == 3);
    assert(tree.root.children[0].name is null);
    assert(tree.root.children[0].children.length == 0);
    assert(tree.root.children[1].name is null);
    assert(tree.root.children[1].children.length == 0);
    child = tree.root.children[2];
    assert(child.parent == tree.root);
    assert(child.name is null);
    assert(child.children.length == 2);
    assert(child.children[0].name is null);
    assert(child.children[1].name is null);
}
