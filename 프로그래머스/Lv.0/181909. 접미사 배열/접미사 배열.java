import java.util.*;
class Solution {
    public String[] solution(String my_string) {
        List<String> list = new ArrayList<>();
        for(int i = 0; i < my_string.length(); i++){
            String curr_string = my_string.substring(i, my_string.length());
            int low = 0;
            int high = list.size();
            while(low < high){
                int mid = (low + high) / 2;
                if(curr_string.compareTo(list.get(mid)) < 0) high = mid;
                else low = mid + 1;
            }
            list.add(low, curr_string);
        }
        
        return list.toArray(new String[0]);
    }
}