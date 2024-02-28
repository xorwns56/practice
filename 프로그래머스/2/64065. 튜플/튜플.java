import java.util.*;
class Solution {
    public int[] solution(String s) {
        String[] sp = s.split("[{}]{1,},?[{}]?");
        sp = Arrays.copyOfRange(sp, 1, sp.length);
        Arrays.sort(sp, (a, b)->{ return a.length() - b.length(); });
        HashSet<Integer> set = new HashSet<>();
        ArrayList<Integer> list = new ArrayList<>();
        for(int i = 0; i < sp.length; i++){
            String[] arr = sp[i].split(",");
            for(int j = 0; j < arr.length; j++){
                int tmp = Integer.parseInt(arr[j]);
                if(!set.contains(tmp)){
                    list.add(tmp);
                    set.add(tmp);
                    break;
                }
            }
        }
        int[] answer = new int[list.size()];
        for(int i = 0; i < list.size(); i++) answer[i] = list.get(i);
        return answer;
    }
}