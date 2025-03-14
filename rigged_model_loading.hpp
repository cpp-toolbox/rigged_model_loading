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

struct VertexBoneData {
    unsigned int indices_of_bones_that_affect_this_vertex[4] = {0};
    float weight_value_of_this_vertex_wrt_bone[4] = {0.0f};

    VertexBoneData() {}

    void add_bone_data(unsigned int BoneID, float Weight);
};

// TODO: move these to draw data?
class IVPNTRigged {
  public:
    IVPNTRigged(std::vector<unsigned int> indices, std::vector<glm::vec3> xyz_positions, std::vector<glm::vec3> normals,
                std::vector<glm::vec2> texture_coordinates, const std::string &texture,
                std::vector<VertexBoneData> bone_data)
        : indices(indices), xyz_positions(xyz_positions), normals(normals), texture_coordinates(texture_coordinates),
          texture_path(texture), bone_data(bone_data) {};
    Transform transform;
    std::vector<unsigned int> indices;
    std::vector<glm::vec3> xyz_positions;
    std::vector<glm::vec3> normals;
    std::vector<glm::vec2> texture_coordinates;
    std::string texture_path;
    std::vector<VertexBoneData> bone_data;
};

// packed version of the above
class IVPNTPRigged {
  public:
    IVPNTPRigged(std::vector<unsigned int> indices, std::vector<glm::vec3> xyz_positions,
                 std::vector<glm::vec3> normals, std::vector<glm::vec2> packed_texture_coordinates,
                 int packed_texture_index, const std::string &texture, std::vector<VertexBoneData> bone_data)
        : indices(indices), xyz_positions(xyz_positions), normals(normals),
          packed_texture_coordinates(packed_texture_coordinates), packed_texture_index(packed_texture_index),
          texture(texture), bone_data(bone_data) {};
    Transform transform;
    int id = GlobalUIDGenerator::get_id();
    std::vector<unsigned int> indices;
    std::vector<glm::vec3> xyz_positions;
    std::vector<glm::vec3> normals;
    std::vector<glm::vec2> packed_texture_coordinates;
    int packed_texture_index;
    std::string texture;
    std::vector<VertexBoneData> bone_data;
};

struct BoneInfo {
    // this transformation puts the local orgin at the start (tail) of the bone and makes the bone sit on an axis so
    // that when transformations are applied it works correctly, note that it is not recursive in any sense, it
    // literally just puts it in the correct position not relative to a parent bone or anything like that

    // another name for this is the inverase bind pose because the bind pose transformation takes a bone and puts it at
    // the origin to be read for application of transformations, and this matrix is the inverse of that

    // A quick thing to jog your memory could be that it "brings the bone joint back to the origin"

    // Note that bones don't really exist, they only exist by the mapping of vertex to bone id's and then knowing the
    // mapping for each bone,

    // Questions for myself: what is a bone tip? it just shows what the rotation and scale is visually, based on a (0,
    // 0, 1) bone tip you can compute where the tip would be based on the scale and so on I think... but then again
    // bones can be larger and have no scaled in them right so then what? it might just be to compute during the auto
    // weighting process maybe check assimp if that data is stored in there if a bone is positioned at position (x, y,
    // z) and then we have a vertex at (x, y + 1, z + 1) then its new position becomes (0, 1, 1) that is it is
    // positioned relative to the bones origin
    glm::mat4 local_space_to_bone_space_in_bind_pose_transformation;
    glm::mat4 local_space_animated_transform_upto_this_bone = glm::mat4(0);

    BoneInfo(const glm::mat4 &lstbst) { local_space_to_bone_space_in_bind_pose_transformation = lstbst; }
};

// NOTE: this is the class you want to use in main
// note that you must make sure all your bone names are unique right now until further improvments
class RecIvpntRiggedCollector {
  public:
    // we make this here so that we always have access to the importer throughout the lifetime of this object
    // without this there will be segfaults with assimp, needs to stay alive
    Assimp::Importer importer;
    const aiScene *scene; // for ease
    std::string directory_to_asset_being_loaded;
    std::vector<IVPNTRigged> ivpntrs;
    int recursion_level_counter = 0;
    void rec_process_nodes(aiNode *node, const aiScene *scene);

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
    std::vector<BoneInfo> bone_unique_idx_to_info;

    // NOTE: ENTRY POINT FUNCTION
    std::vector<IVPNTRigged> parse_model_into_ivpntrs(const std::string &model_path);
    std::vector<VertexBoneData> process_mesh_vertices_bone_data(aiMesh *mesh);
    IVPNTRigged process_mesh_ivpntrs(aiMesh *mesh, const aiScene *scene);

    unsigned int curr_animation_index_rec = 0;
    // the armature node name is whatever the node for the armature is called in blender
    // NOTE: this is populated during the initial parse model function, it is used during rec_update_animation_matrices
    //     - we recursively iterate over the hierarchical structure of the assimp import
    //     - if the node is an armature, then we parse
    std::unordered_map<std::string, std::unordered_map<std::string, int>>
        armature_node_name_to_animation_name_to_assimp_animation_index;

    // this is used for for eventually binding into uniforms with all the matrices, then in the shader
    // we also have a vertex attribute for each vertex which specifies the id of which matrices to use...
    void set_bone_transforms(float time_in_seconds, std::vector<glm::mat4> &transforms_to_be_set,
                             std::string requested_animation);
    void update_animation_matrices(float animation_time_ticks, std::string requested_animation);
    void rec_update_animation_matrices(float animation_time_ticks, glm::mat4 parent_transform, aiNode *node,
                                       const aiScene *scene, int rec_depth, std::string requested_animation);
};

const aiNodeAnim *find_node_anim(const aiAnimation *pAnimation, const std::string &NodeName);
std::string get_full_node_path(const aiNode *node);
void print_all_animations(const aiScene *scene);
void print_ai_animation(const aiAnimation *anim);
unsigned int find_animation_index_by_name(const aiScene *scene, const std::string &animationName);

bool is_armature_node(const aiNode *node);
std::unordered_map<std::string, std::unordered_map<std::string, int>>
build_armature_name_to_animation_name_to_assimp_animation_index_map(const aiScene *scene);

unsigned int find_idx_of_scaling_key_for_given_time(float animation_time_ticks, const aiNodeAnim *node_anim);
unsigned int find_idx_of_rotation_key_for_given_time(float animation_time_ticks, const aiNodeAnim *node_anim);
unsigned int find_idx_of_translation_key_for_given_time(float animation_time_ticks, const aiNodeAnim *node_anim);
void calc_interpolated_scaling(aiVector3D &out, float animation_time_ticks, const aiNodeAnim *node_anim);
void calc_interpolated_rotation(aiQuaternion &out, float animation_time_ticks, const aiNodeAnim *node_anim);
void calc_interpolated_translation(aiVector3D &out, float animation_time_ticks, const aiNodeAnim *node_anim);

} // namespace rigged_model_loading

#endif // RIGGED_MODEL_LOADING_HPP
