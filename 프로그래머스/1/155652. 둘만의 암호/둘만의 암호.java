class Solution {
    public String solution(String s, String skip, int index) {
        char[] chars = s.toCharArray();
        boolean[] isSkip = new boolean[26];
        for(char c : skip.toCharArray()) isSkip[c - 'a'] = true;
        for(int i = 0; i < chars.length; i++){
            int count = 0;
            int curr_idx = chars[i] - 'a';
            while(count < index){
                curr_idx = (curr_idx + 1) % 26;
                if(!isSkip[curr_idx]) count++;
            }
            chars[i] = (char)('a' + curr_idx);
        }
        return String.valueOf(chars);
    }
}