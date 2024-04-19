import java.util.*;
public class File implements Comparable<File>{
    String head;
    String number;
    String tail;
    @Override
    public String toString(){
        return head + number + tail;
    }
    @Override
    public int compareTo(File file){
        int comp = head.toUpperCase().compareTo(file.head.toUpperCase());
        if(comp != 0) return comp;
        return Integer.parseInt(number) - Integer.parseInt(file.number);
    }
}
class Solution {
    public String[] solution(String[] files) {
        File[] arr = new File[files.length];
        for(int i = 0; i < files.length; i++){
            int idx = 0;
            arr[i] = new File();
            while(idx < files[i].length()){
                char c = files[i].charAt(idx);
                if('0' <= c && c <= '9') break;
                idx++;
            }
            arr[i].head = files[i].substring(0, idx);
            int tmp = idx;
            while(idx < files[i].length()){
                char c = files[i].charAt(idx);
                if(!('0' <= c && c <= '9')) break;
                idx++;
            }
            arr[i].number = files[i].substring(tmp, idx);
            arr[i].tail = files[i].substring(idx, files[i].length());
        }
        Arrays.sort(arr);
        String[] answer = new String[files.length];
        for(int i = 0; i < files.length; i++) answer[i] = arr[i].toString();
        return answer;
    }
}