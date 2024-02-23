import java.util.*;
class Solution {
    public int[][] solution(int[][] data, String ext, int val_ext, String sort_by) {
        HashMap<String, Integer> map = new HashMap<>();
        map.put("code", 0);
        map.put("date", 1);
        map.put("maximum", 2);
        map.put("remain", 3);
        List<int[]> filtered_list = new ArrayList<>();
        int ext_idx = map.get(ext);
        int sort_by_idx = map.get(sort_by);
        for(int i = 0; i < data.length; i++){
            if(data[i][ext_idx] < val_ext) filtered_list.add(data[i]);
        }
        Collections.sort(filtered_list, (a, b)->{
            return a[sort_by_idx] - b[sort_by_idx];
        });
        return filtered_list.toArray(new int[0][]);
    }
}