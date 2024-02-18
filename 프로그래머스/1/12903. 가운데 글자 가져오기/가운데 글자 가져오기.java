class Solution {
    public String solution(String s) {
        int start = s.length() / 2;
        int end = start + 1;
        if((s.length() & 1) == 0) start--;
        return s.substring(start, end);
    }
}