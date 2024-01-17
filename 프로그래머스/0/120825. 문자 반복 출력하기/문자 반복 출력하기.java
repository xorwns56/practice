class Solution {
    public String solution(String my_string, int n) {
        char[] chars = new char[my_string.length() * n];
        for(int i = 0; i < chars.length; i++) chars[i] = my_string.charAt(i / n);
        return String.valueOf(chars);
    }
}