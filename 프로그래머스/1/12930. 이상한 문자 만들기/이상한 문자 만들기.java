class Solution {
    public String solution(String s) {
        char[] chars = s.toCharArray();
        int n = 0;
        for(int i = 0; i < chars.length; i++){
            if(chars[i] != ' '){
                if((n & 1) == 0 && 'a' <= chars[i] && chars[i] <= 'z') chars[i] += 'A' - 'a';
                else if((n & 1) == 1 && 'A' <= chars[i] && chars[i] <= 'Z') chars[i] += 'a' - 'A';
                n++;
            }else n = 0;
        }
        return String.valueOf(chars);
    }
}