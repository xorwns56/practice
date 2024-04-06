import java.util.*;
class Solution {
    public int[] solution(String[] operations) {
        List<Integer> list = new ArrayList<>();
        for(int i = 0; i < operations.length; i++){
            String[] sp = operations[i].split("\\s");
            if(sp[0].equals("I")){
                int value = Integer.parseInt(sp[1]);
                int low = 0;
                int high = list.size();
                while(low < high){
                    int mid = (low + high) / 2;
                    if(list.get(mid) < value) low = mid + 1;
                    else high = mid;
                }
                list.add(low, value);
            }else if(list.size() > 0){
                if(sp[1].charAt(0) == '-') list.remove(0);
                else list.remove(list.size() - 1);
            }
        }
        int[] answer = list.size() == 0 ? new int[] { 0, 0 } : new int[] { list.get(list.size() - 1), list.get(0) };
        return answer;
    }
}