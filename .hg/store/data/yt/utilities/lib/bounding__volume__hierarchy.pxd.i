        �  �      J<���������P>��Yf�Oܖ�%{FJҎc            x��S�n�0��n�n��n�ةS�CѥE�m�dH��}iǩ[�U��H��=�����Tkg!���=6=� ��T�~��=h%��[[���XNW�'�">��/�����y�Ⱦ�ޯC�<�H.%�	�HW2�I��DaQ�^�0�3�����M�(�IC�E�ZM��Ӎ��GwJh4XF���$Y6�����:�P:�ƹW�L{D�	����o�>=sU�Tp@6w%�����#�[d�t�X�Í�EO�V��˔��yH�����+x9�U�B�#}�����)#CD�y��,ٹ8���a^7}+r(��:���"�B��~�M��F�/擏�|~���)�u!4��f�ѱ#�-}�u��2n��]�
yE�1�6%x4r������:�>P��3����ġ%����&{���4    �     �        J?    ����,��;���J�U�K�e,P�            x�c``ng``����P����e敘�ė($V耸i9���ₜ�M�����.����4������Ԣ�����Ԝ4���J-���J�je ~
t�b�OJM���AI�K�n}|QjriQqfYj<.��(8�y�姤jMHI�0.�T����M�s 7Vd�    G     `       JA   ������*���gjT�=EJj��Ċ��            x�c``vd``>���Я �)�i
yz�y%f&�%
��E%�%��yũ9i:�rI��y("�y):\
����
��D�@qANf��B^~zf �9*�    �     (  !     JE   �����*������s_���	�0�h�              �  �       np.float64_t d0, d1, d2
    �     n  d     JH   �����ߎ�.!����a�$wR�M            x�c``��������Х �)�i
Na~�)�Z
�E�ɥEře��I��9)ũ9i:
yz�y%f&�%
I��y("�y)�
y��9\p��3S�K+JM���+H�4 QK/�    =     *  �     JI   ����Y>���B�Q�\-$g\mj              6  6       cdef np.int64_t leaf_size
    g       =     J�   ������(�t��u~�t�}�)8t�              6  {        s     N  �     M�   ����ܷhϼ���T�';=���M            x�c``6c c� HNIMS�+�K��O,13�/���Q��2�U(K-*�LN-�BV����(3/C�Qi��9)�)�%�\ cQ$7    �       �     M�   ���� A�+LwU�$[�����,�[|              �  �        �     �  F     M�   ����.j_�z��Z��j=..�i	            x�}�1�0E#��L�2!��8�B� Kq%n%.�5�)e�Dk�O�����Yͺ�����@BXu��ƀ�Qۜ)>*8✊�u7�9`���^-��\�b�>$+����X��n�{�Лފ]8(~9ā�2��Zm�V���T���_N�����;f�    l     )  c   	  M�   	������r����Uij������(��              Q  Q       np.int64_t near_boundary
    �     3  �   
  M�   
�����۞�G� ��%�K�v��              <  <   '    cdef np.int64_t num_field_per_elem
    �     G  �     N6   ����O6�z�,�Ra� �>���            x�c``fe``6` �Rf�P ���4�������3��--��Լ���̔b.�
'��
-������b. =4�         J  <     N7   �����H�O[�$��zb#gUQH �               u   u   2from yt.utilities.lib.primitives cimport Triangle
  y          Y     �  	|     N:   ������[����V�BjA��Ǧ��            x���Mj�0�E�J��&��� G�6��cWk�5��er�J
�[�E�V0 �<=͇$��C��g��|G�2v�uo5���%�����\gZ�怐��IjVX�uO���R0�U�1����\��2���R���9��U�����.]9,�Rc�r"�F�] �#x
ȗ�mbPS����{*ӱ��6}'�z[ۄ���wA��;�W+:*�Cd<!^�%����~h���,��'O2Ei4���ɕu�p%��    V     J  	z     N;   ����v����%�=y�b5���5�x            x�c``��������AJ�SR���3S�
�2s3K2�R���"���J�*�@/3���$�D!�47�� �N    �    @  '     N<   �����>�p.�OG�j%AF��_LI            x����J�0�SA�>��]�?��(^�ʁE���N����6Ц%I�|#����YL��U�X��9'�/|IBo��5#�1l�cQ�x��LIXp�D*V�y:ql�x��U]��3yx\ouo7C��s�p�����s�բ��vē{s���v�#�`�'UF��~/F�_��%�(��~�6B�#d�ͻ:��h���s}��%Qg�,',��<8�Zkr:@b����K2�kǇ:���M��U[����β)X��E���f O$d����J'���6}g�� o@(�en��Z�j�)�E�2{8bIU����x��qo��    �     &  A     N=   ����|�kA���fe�VQ�T_�����              �  �       int hex20_faces[6][8]
         K  
�     NG   ����(Hۈj��G�g�:��l܆               u   �   3from yt.utilities.lib.primitives cimport BBox, Ray
  �  M        Q     5  
�     N�   ������8�d����Q�@["���               �   �   )cdef extern from "mesh_triangulation.h":
    �     s  
�     SW   �����G�pB�fQ��O{JG            x�c`  �O��-�/*QH�,����b``|	�R
@��W�P�Zbh����ZmmT�9��U�T%���)��g�(ħ����姤�'%�Wh���(8�y���@:\ ~a!�    �     g  �     SX   ����F<x�3 BqS��~l8���#�3            x�c`�g c��
@�����P����_�Z_Z_�
K�3R�5�Ss�t���������3���h+�X��Ԣ����b"�f�!k��K��T��O��� �*+)