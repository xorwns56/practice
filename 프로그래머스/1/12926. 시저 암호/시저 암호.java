class Solution {
    public String solution(String s, int n) {
        char[] chars = s.toCharArray();
        for(int i = 0; i < chars.length; i++){
            if(chars[i] != ' '){
                if('a' <= chars[i] && chars[i] <= 'z') chars[i] = (char)('a' + (chars[i] + n - 'a') % 26);
                else if('A' <= chars[i] && chars[i] <= 'Z') chars[i] = (char)('A' + (chars[i] + n - 'A') % 26);
            }
        }
        return String.valueOf(chars);
    }
}