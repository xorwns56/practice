class Solution {
    public int solution(int[] bandage, int health, int[][] attacks) {
        int last_attack_frame = attacks[attacks.length - 1][0];
        int attack_idx = 0;
        int curr_health = health;
        int heal_count = 0;
        for(int frame = 0; frame <= last_attack_frame; frame++){
            if(attacks[attack_idx][0] == frame){
                heal_count = 0;
                curr_health -= attacks[attack_idx][1];
                if(curr_health <= 0) return -1;
                attack_idx++;
            }else{
                heal_count++;
                curr_health += bandage[1];
                if(heal_count == bandage[0]){
                    curr_health += bandage[2];
                    heal_count = 0;
                }
                curr_health = Math.min(curr_health, health);
            }
        }
        return curr_health;
    }
}