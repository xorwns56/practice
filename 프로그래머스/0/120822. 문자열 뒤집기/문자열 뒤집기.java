class Solution {
    public String solution(String my_string) {
        char[] chars = my_string.toCharArray();
        for(int i = 0; i < chars.length / 2; i++){
            char tmp = chars[i];
            chars[i] = chars[chars.length - 1 - i];
            chars[chars.length - 1 - i] = tmp;
        }
        return String.valueOf(chars);
    }
}