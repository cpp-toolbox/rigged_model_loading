#ifndef RIGGED_MODEL_LOADING_HPP
#define RIGGED_MODEL_LOADING_HPP

#include <glm/glm.hpp>
#include <string>
#include <utility>
#include <vector>
#include <assimp/scene.h>
#include <assimp/Importer.hpp>
#include <functional>

#include "sbpt_generated_includes.hpp"

namespace rigged_model_loading {

// NOTE: this is the class you want to use in main
// note that you must make sure all your bone names are unique right now until further improvments
class RecIvpntRiggedCollector {
  public:
    RecIvpntRiggedCollector(UniqueIDGenerator &ivpntr_id_generator) : ivpntr_id_generator(ivpntr_id_generator) {};
    // we make this here so that we always have access to the importer throughout the lifetime of this object
    // without this there will be segfaults with assimp, needs to stay alive
    Assimp::Importer importer;
    const aiScene *scene;                   // for ease
    UniqueIDGenerator &ivpntr_id_generator; // used to give each ivpntr an id for rendering
    std::string directory_to_asset_being_loaded;
    std::vector<draw_info::IVPNTRigged> ivpntrs;
    int recursion_level_counter = 0;
    int no_anim_sentinel = -1;
    std::string current_animation_name;
    bool animation_is_complete = false;
    double current_animation_time = 0;
    void rec_process_nodes(aiNode *node, const aiScene *scene);
    glm::mat4 get_the_transform_to_attach_an_object_to_a_bone(std::string bone_name,
                                                              Transform &bone_origin_attachment_offset);

    // we work with the transforms of the bones of a mesh, but we never look at the root bone/node's transform
    // if it does have a transform, then in that case it means that the model is already placed somewhere in the world,
    // therefore if we only used the transforms of the non-root nodes it would never get placed where it has to be,
    // therefore since this transform shows how to go from the root node to the origin, then the inverse takes us from
    // the origin to the root node.
    glm::mat4 inverse_root_node_transform;
    int get_next_bone_id(const aiBone *pBone);
    // NOTE: this is populated during the initial parse_model function
    // bone names are parsed from the file directly
    // TODO: use a id generator in the future ?
    std::unordered_map<std::string, int> bone_name_to_unique_index; // one for each bone.
    std::vector<draw_info::BoneInfo> bone_unique_idx_to_info;

    // NOTE: ENTRY POINT FUNCTION
    std::vector<draw_info::IVPNTRigged> parse_model_into_ivpntrs(const std::string &model_path);
    std::vector<draw_info::VertexBoneData> process_mesh_vertices_bone_data(aiMesh *mesh);
    draw_info::IVPNTRigged process_mesh_ivpntrs(aiMesh *mesh, const aiScene *scene);

    unsigned int curr_animation_index_rec = 0;
    // the armature node name is whatever the node for the armature is called in blender
    // NOTE: this is populated during the initial parse model function, it is used during rec_update_animation_matrices
    //     - we recursively iterate over the hierarchical structure of the assimp import
    //     - if the node is an armature, then we parse
    std::unordered_map<std::string, std::unordered_map<std::string, int>>
        armature_node_name_to_animation_name_to_assimp_animation_index;

    // note that an animation might have multiple animation indices, because you
    // get one for each armature involved in an animation, for our purposes, we simply use
    // any such one of those armature indices, we run under the assumption that the metadata
    // about each one is the same (which might be the cause of some errors)
    std::unordered_map<std::string, int> animation_name_to_assimp_animation_index;

    std::unordered_map<std::string, std::unordered_map<std::string, int>>
    build_armature_name_to_animation_name_to_assimp_animation_index_map(const aiScene *scene);

    // this is used for for eventually binding into uniforms with all the matrices, then in the shader
    // we also have a vertex attribute for each vertex which specifies the id of which matrices to use...
    void set_bone_transforms(float delta_time, std::vector<glm::mat4> &transforms_to_be_set,
                             std::string requested_animation, bool loop = false, bool restart = false,
                             bool hold_last_frame = false);
    void update_animation_matrices(float animation_time_ticks, std::string requested_animation);
    void rec_update_animation_matrices(float animation_time_ticks, glm::mat4 parent_transform, aiNode *node,
                                       const aiScene *scene, int rec_depth, std::string requested_animation);
};

const aiNodeAnim *find_node_anim(const aiAnimation *pAnimation, const std::string &NodeName);
std::string get_full_node_path(const aiNode *node);
void print_all_animations(const aiScene *scene);
void print_ai_animation_short(const aiAnimation *anim);
void print_ai_animation(const aiAnimation *anim);
void print_ai_node_anim(const aiNodeAnim *anim);
unsigned int find_animation_index_by_name(const aiScene *scene, const std::string &animationName);

bool is_armature_node(const aiNode *node);

unsigned int find_idx_of_scaling_key_for_given_time(float animation_time_ticks, const aiNodeAnim *node_anim);
unsigned int find_idx_of_rotation_key_for_given_time(float animation_time_ticks, const aiNodeAnim *node_anim);
unsigned int find_idx_of_translation_key_for_given_time(float animation_time_ticks, const aiNodeAnim *node_anim);
void calc_interpolated_scaling(aiVector3D &out, float animation_time_ticks, const aiNodeAnim *node_anim);
void calc_interpolated_rotation(aiQuaternion &out, float animation_time_ticks, const aiNodeAnim *node_anim);
void calc_interpolated_translation(aiVector3D &out, float animation_time_ticks, const aiNodeAnim *node_anim);

} // namespace rigged_model_loading

#endif // RIGGED_MODEL_LOADING_HPP
