class Solution {
    public String solution(String my_string, int s, int e) {
        char[] chars = my_string.toCharArray();
        int s_idx = s;
        int e_idx = e;
        while(s_idx < e_idx){
            char tmp = chars[s_idx];
            chars[s_idx] = chars[e_idx];
            chars[e_idx] = tmp;
            s_idx++;
            e_idx--;
        }
        return String.valueOf(chars);
    }
}