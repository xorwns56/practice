class Solution {
    public String solution(String s) {
        char[] chars = s.toCharArray();
        boolean first = true;
        for(int i = 0; i < chars.length; i++){
            if(chars[i] != ' '){
                if(first && 'a' <= chars[i] && chars[i] <= 'z') chars[i] = (char)(chars[i] - 'a' + 'A');
                else if(!first && 'A' <= chars[i] && chars[i] <= 'Z') chars[i] = (char)(chars[i] - 'A' + 'a');
                first = false;
            }else first = true;
        }
        return String.valueOf(chars);
    }
}