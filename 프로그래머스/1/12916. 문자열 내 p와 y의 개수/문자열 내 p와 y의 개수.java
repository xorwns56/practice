class Solution {
    boolean solution(String s) {
        char[] chars = s.toLowerCase().toCharArray();
        int count = 0;
        for(int i = 0; i < chars.length; i++){
            if(chars[i] == 'p') count++;
            else if(chars[i] == 'y') count--;
        }
        return count == 0;
    }
}