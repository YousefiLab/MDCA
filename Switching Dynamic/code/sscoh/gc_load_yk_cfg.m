function []= gc_load_yk_cfg(ch_num)

str_dict = 'C:\Users\REZA_SF\Desktop\validation\code\GC_empirical_with_slepian_anesthesia\';

switch ch_num
    case 32
%         str_name_1 = sprintf('Anethesia_Yk_chNum%d_fr11_sr256_winSec8_segNum8_TRIAL', ch_num);
%         str_name_2 = sprintf('Anethesia_cfg_Init_chNum%d_fr11_sr256_winSec8_segNum8_TRIAL', ch_num);
        
        str_name_1 = sprintf('Anethesia_Yk_chNum%d_fr13_sr256_winSec32_segNum32_TRIAL', ch_num);
        str_name_2 = sprintf('Anethesia_cfg_Init_chNum%d_fr13_sr256_winSec32_segNum32_TRIAL', ch_num);
        
    case 20
        str_name_1 = sprintf('Anethesia_Yk_chNum%d_fr13_sr256_winSec20_segNum20_TRIAL', ch_num);
        str_name_2 = sprintf('Anethesia_cfg_Init_chNum%d_fr13_sr256_winSec20_segNum20_TRIAL', ch_num);
        
    case 10
%         str_name_1 = sprintf('Anethesia_Yk_chNum%d_fr12_sr256_winSec20_segNum10_TRIAL', ch_num);
%         str_name_2 = sprintf('Anethesia_cfg_Init_chNum%d_fr12_sr256_winSec20_segNum10_TRIAL', ch_num);
        
        str_name_1 = sprintf('Anethesia_Yk_chNum%d_fr13_sr256_winSec20_segNum20_TRIAL', ch_num);
        str_name_2 = sprintf('Anethesia_cfg_Init_chNum%d_fr13_sr256_winSec20_segNum20_TRIAL', ch_num);      
     
    case 8
        str_name_1 = sprintf('Anethesia_cfg_Init_chNum%d_fr13_sr256_winSec16_segNum16_TRIAL', ch_num);
        str_name_2 = sprintf('Anethesia_Yk_chNum%d_fr13_sr256_winSec16_segNum16_TRIAL', ch_num);
             
    case 5
        str_name_1 = sprintf('Anethesia_Yk_chNum%d_fr13_sr256_winSec8_segNum8_TRIAL', ch_num);
        str_name_2 = sprintf('Anethesia_cfg_Init_chNum%d_fr13_sr256_winSec8_segNum8_TRIAL', ch_num);
        
    case 3
        str_name_1 = sprintf('Anethesia_Yk_chNum%d_fr13_sr256_winSec15_segNum15_TRIAL', ch_num);
        str_name_2 = sprintf('Anethesia_cfg_Init_chNum%d_fr13_sr256_winSec15_segNum15_TRIAL', ch_num);
        
%         str_name_1 = sprintf('Anethesia_Yk_chNum%d_fr13_sr256_winSec10_segNum10_TRIAL', ch_sim);
%         str_name_2 = sprintf('Anethesia_cfg_Init_chNum%d_fr13_sr256_winSec10_segNum10_TRIAL', ch_sim);
end

str_load_1 = sprintf('%s%s', str_dict, str_name_1);
str_load_2 = sprintf('%s%s', str_dict, str_name_2);


load(str_load_1)
load(str_load_2)

end