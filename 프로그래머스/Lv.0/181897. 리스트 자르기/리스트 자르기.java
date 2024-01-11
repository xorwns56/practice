import java.util.*;
class Solution {
    public List<Integer> solution(int n, int[] slicer, int[] num_list) {
        List<Integer> list = new ArrayList<>();
        int start = slicer[0];
        int end = slicer[1];
        if(n == 1) start = 0;
        else if(n == 2) end = num_list.length - 1;
        for(int i = start; i <= end; i += n == 4 ? slicer[2] : 1) list.add(num_list[i]);
        return list;
    }
}