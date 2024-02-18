class Solution {
    public boolean solution(String s) {
        char[] chars = s.toCharArray();
        if(chars.length != 4 && chars.length != 6) return false;
        for(int i = 0; i < chars.length; i++){
            if(!('0' <= chars[i] && chars[i] <= '9')) return false;
        }
        return true;
    }
}