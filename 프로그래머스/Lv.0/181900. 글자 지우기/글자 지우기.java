class Solution {
    public String solution(String my_string, int[] indices) {
        char[] chars = my_string.toCharArray();
        for(int i = 0; i < indices.length; i++) chars[indices[i]] = ' ';
        return String.valueOf(chars).replaceAll(" ", "");
    }
}