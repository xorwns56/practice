class Solution {
    public String solution(String X, String Y) {
        int[][] count = new int[2][10];
        for(char c : X.toCharArray()) count[0][c - '0']++;
        for(char c : Y.toCharArray()) count[1][c - '0']++;
        StringBuilder sb = new StringBuilder();
        for(int i = 9; i >= 0; i--){
            while(count[0][i] > 0 && count[1][i] > 0){
                boolean escape = i == 0 && sb.length() == 0;
                sb.append(i);
                if(escape) break;
                count[0][i]--;
                count[1][i]--;
            }
        }
        return sb.length() == 0 ? "-1" : sb.toString();
    }
}