import java.util.*;
class Solution {
    public String[] solution(String myString) {
        List<String> list = new ArrayList<>();
        String tmp = "";
        for(int i = 0; i < myString.length(); i++){
            char c = myString.charAt(i);
            if(c == 'x' && !tmp.isEmpty()){
                insert(list, tmp);
                tmp = "";
            }else if(c != 'x') tmp += c;
        }
        if(!tmp.isEmpty()) insert(list, tmp);
        return list.toArray(new String[0]);
    }
    
    public void insert(List<String> list, String str){
        int low = 0;
        int high = list.size();
        while(low < high){
            int mid = (low + high) / 2;
            if(str.compareTo(list.get(mid)) < 0) high = mid;
            else low = mid + 1;
        } 
        list.add(low, str);
    }
}