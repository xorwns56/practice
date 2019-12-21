/*
We run a preorder depth first search on the root of a binary tree.
At each node in this traversal, we output D dashes (where D is the depth of this node), then we output the value of this node.  (If the depth of a node is D, the depth of its immediate child is D+1.  The depth of the root node is 0.)
If a node has only one child, that child is guaranteed to be the left child.
Given the output S of this traversal, recover the tree and return its root.

Example 1:
Input: "1-2--3--4-5--6--7"
Output: [1,2,5,3,4,6,7]

Example 2:
Input: "1-2--3---4-5--6---7"
Output: [1,2,5,3,null,6,null,4,null,7]
 

Example 3:
Input: "1-401--349---90--88"
Output: [1,401,null,349,88,90]
 
Note:
The number of nodes in the original tree is between 1 and 1000.
Each node will have a value between 1 and 10^9.

 * Definition for a binary tree node.
 * public class TreeNode {
 *     int val;
 *     TreeNode left;
 *     TreeNode right;
 *     TreeNode(int x) { val = x; }
 * }
*/
class Solution {
    public TreeNode recoverFromPreorder(String S) {
        char[] c = S.toCharArray();
        int depth = 0;
        String numStr = "";
        List<TreeNode> treeNodeList = new ArrayList<TreeNode>();
        for(int i=0;i<c.length;i++){
            if(c[i]=='-') depth++;
            else{
                numStr += c[i];
                if(i+1>=c.length||c[i+1]=='-'){
                    TreeNode treeNode = new TreeNode(Integer.parseInt(numStr));
                    if(depth>=treeNodeList.size()) treeNodeList.add(treeNode);
                    else{
                        treeNodeList.remove(depth);
                        treeNodeList.add(depth, treeNode);
                    }
                    
                    if(depth>0){
                        if(treeNodeList.get(depth-1).left==null) treeNodeList.get(depth-1).left = treeNode;
                        else treeNodeList.get(depth-1).right = treeNode;
                    }
                    
                    depth = 0;
                    numStr = "";
                }
            }
        }
        
        return treeNodeList.get(0);
    }
}
