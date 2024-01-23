class Solution {
    public String[] solution(String[] quiz) {
        String[] answer = new String[quiz.length];
        for(int i = 0; i < quiz.length; i++){
            String[] sp = quiz[i].replace("- ", "-").replace("--", "").replace("+ ", "").replace("= ", "").split("\\s");
            int sum = 0;
            for(int j = 0; j < sp.length; j++){
                if(j == sp.length - 1) answer[i] = sum == Integer.parseInt(sp[j]) ? "O" : "X";
                else sum += Integer.parseInt(sp[j]);
            }
        }
        return answer;
    }
}