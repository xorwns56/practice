class Solution {
    public int[] solution(String my_string) {
        int[] answer = new int[52];
        for(char c : my_string.toCharArray()) answer['A' <= c && c <= 'Z' ? (c - 'A') : (c + 26 - 'a')]++;
        return answer;
    }
}