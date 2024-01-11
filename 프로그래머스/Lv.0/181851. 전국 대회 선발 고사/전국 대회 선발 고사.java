import java.util.*;
class Solution {
    public int solution(int[] rank, boolean[] attendance) {
        List<Integer> idx_list = new ArrayList<>();
        for(int i = 0; i < rank.length; i++){
            if(!attendance[i]) continue;
            int low = 0;
            int high = idx_list.size();
            while(low < high){
                int mid = (low + high) / 2;
                if(rank[i] < rank[idx_list.get(mid)]) high = mid;
                else low = mid + 1;
            }
            idx_list.add(low, i);
        }
        return 10000 * idx_list.get(0) + 100 * idx_list.get(1) + idx_list.get(2);
    }
}