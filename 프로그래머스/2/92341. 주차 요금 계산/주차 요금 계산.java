import java.util.*;
class Solution {
    public int[] solution(int[] fees, String[] records) {
        TreeMap<String, String> map = new TreeMap<>();
        for(int i = 0; i < records.length; i++){
            String[] sp = records[i].split("\\s");
            String prev = map.get(sp[1]);
            map.put(sp[1], prev != null ? prev + " " + sp[0] : sp[0]);
        }
        Iterator<String> iter = map.values().iterator();
        int[] answer = new int[map.size()];
        for(int i = 0; i < answer.length; i++){
            String[] sp = iter.next().split("\\s");
            if(sp.length % 2 == 1){
                sp = Arrays.copyOfRange(sp, 0, sp.length + 1);
                sp[sp.length - 1] = "23:59";
            }
            int total = 0;
            for(int j = 0; j < sp.length; j += 2){
                String[] in = sp[j].split(":");
                String[] out = sp[j + 1].split(":");
                total += Integer.parseInt(out[0]) * 60 + Integer.parseInt(out[1]) - (Integer.parseInt(in[0]) * 60 + Integer.parseInt(in[1]));
            }
            answer[i] = fees[1];
            if(total > fees[0]) answer[i] += (total - fees[0] + fees[2] - 1) / fees[2] * fees[3];
        }
        return answer;
    }
}