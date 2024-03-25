import java.util.*;
class Solution {
    public int[][] solution(int n) {
        List<int[]> list = new LinkedList<>();
        hanoi(list, n, 1, 2, 3);
        return list.toArray(new int[list.size()][]);
    }
    public void hanoi(List<int[]> list, int n, int from, int by, int to){
        if(n == 0) return;
        hanoi(list, n - 1, from, to, by);
        list.add(new int[] { from, to });
        hanoi(list, n - 1, by, from, to);
    }
}