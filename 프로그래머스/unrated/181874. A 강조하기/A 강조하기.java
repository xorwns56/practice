class Solution {
    public String solution(String myString) {
        char[] chars = myString.toCharArray();
        for(int i = 0; i < chars.length; i++){
            if(chars[i] == 'a') chars[i] += 'A' - 'a';
            else if('B' <= chars[i] && chars[i] <= 'Z') chars[i] += 'a' - 'A';
        }
        return String.valueOf(chars);
    }
}