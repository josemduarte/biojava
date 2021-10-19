package org.biojava.nbio.structure.jaligner;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.StructureAlignment;
import org.biojava.nbio.structure.align.StructureAlignmentFactory;
import org.biojava.nbio.structure.align.fatcat.FatCatRigid;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.jama.Matrix;
import org.glassfish.jersey.media.multipart.FormDataMultiPart;
import org.glassfish.jersey.media.multipart.MultiPartFeature;
import org.junit.Test;

import javax.vecmath.Matrix4d;
import javax.vecmath.Vector3d;
import javax.ws.rs.client.Client;
import javax.ws.rs.client.ClientBuilder;
import javax.ws.rs.client.Entity;
import javax.ws.rs.core.MediaType;
import javax.ws.rs.core.Response;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class TestAlignmentFromTransform {

    private static final String BASE_URL = "https://alignment-west.rcsb.org";
    private static final String SUBMISSION_ENDPOINT = "/api/v1-beta/structures/submit";
    private static final String RESULTS_ENDPOINT = "/api/v1-beta/structures/results?uuid=";

    private static final Matrix4d ID_MATRIX = new Matrix4d(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1);

    //-F query='{"mode":"pairwise","method":{"name":"fatcat-rigid"},"structures":[{"entry_id":"1FSL","asym_id":"A"},{"entry_id":"4HHB","asym_id":"A"}]}'

    @Test
    public void testAlignerAgainstRcsbAlignmentService() throws Exception {

        String pdbId1 = "1FSL";
        String asymId1 = "A";
        String pdbId2 = "4HHB";
        String asymId2 = "A";
        StrucAlignment strucAli = getAlignmentFromRcsbAlignmentService(pdbId1, asymId1, pdbId2, asymId2);

        Structure s1 = StructureIO.getStructure(pdbId1);
        Structure s2 = StructureIO.getStructure(pdbId2);

        Chain c1 = s1.getPolyChain(asymId1);
        Chain c2 = s2.getPolyChain(asymId2);

        StrucAlignment strucAliSelf = computeFatcatAlignment(c1, c2);

        Calc.transform(c1, strucAli.matrices.get(0));
        Calc.transform(c2, strucAli.matrices.get(1));

        AtomGroupSequence atoms1 = new AtomGroupSequence(c1.getAtomGroups());
        AtomGroupSequence atoms2 = new AtomGroupSequence(c2.getAtomGroups());
        long start = System.currentTimeMillis();
        Alignment ali = NeedlemanWunschGotoh.align(atoms1, atoms2, 6, 1);
        long end = System.currentTimeMillis();
        System.out.println("Alignment calculated in " + (end-start) + "ms");
        System.out.println(ali.getSequence1());
        System.out.println(ali.getMarkupLine());
        System.out.println(ali.getSequence2());
    }

    private static class StrucAlignment {
        List<Matrix4d> matrices;
        List<String> seqs;
        private StrucAlignment(List<Matrix4d> matrices, List<String> seqs) {
            this.matrices = matrices;
            this.seqs = seqs;
        }
    }

    private StrucAlignment getAlignmentFromRcsbAlignmentService(String entryId1, String asymId1, String entryId2, String asymId2) throws IOException, InterruptedException {
        String query = String.format(
                "{" +
                    "\"context\":{" +
                    "\"mode\":\"pairwise\"," +
                    "\"method\":{\"name\":\"fatcat-rigid\"}," +
                    "\"structures\":[{\"entry_id\":\"%s\",\"asym_id\":\"%s\"},{\"entry_id\":\"%s\",\"asym_id\":\"%s\"}]" +
                "}}", entryId1, asymId1, entryId2, asymId2);

        FormDataMultiPart multiPart = new FormDataMultiPart();
        multiPart.setMediaType(MediaType.MULTIPART_FORM_DATA_TYPE);
        multiPart.field("query", query);

        Client client = ClientBuilder.newBuilder()
                .register(MultiPartFeature.class)
                .build();
        Response response = client.target(BASE_URL + SUBMISSION_ENDPOINT)
                .request(MediaType.TEXT_PLAIN)
                .header("Content-Type", "multipart/form-data")
                .post(Entity.entity(multiPart, multiPart.getMediaType()));
        if (response.getStatus() !=200)
            throw new IOException("Bad response status: " + response.getStatus());
        String token = response.readEntity(String.class);

        JsonNode node = null;
        long start = System.currentTimeMillis();
        while (System.currentTimeMillis() - start < 30000) {
            response = client.target(BASE_URL + RESULTS_ENDPOINT + token)
                    .request(MediaType.APPLICATION_JSON).get();
            String jsonResp = response.readEntity(String.class);
            ObjectMapper mapper = new ObjectMapper();
            node = mapper.readValue(jsonResp, JsonNode.class);

            String status = node.get("info").get("status").asText();
            if (status.equals("COMPLETE")) {
                break;
            } else {
                Thread.sleep(500);
            }
            node = null;
        }
        if (node == null) throw new IOException("Timed out");
        //System.out.println(node.get("results").toString());
        List<Matrix4d> matrices = new ArrayList<>();
        List<String> seqs = new ArrayList<>();
        for (JsonNode oneResult : node.get("results")) {
            for (JsonNode oneBlock : oneResult.get("blocks")) {
                for (JsonNode t : oneBlock.get("transformations")) {
                    Matrix4d m = new Matrix4d();
                    int i = 0, j = 0;
                    for (JsonNode v : t) {
                        m.setElement(i++, j, v.asDouble());
                        if (i == 4) {
                            j++;
                            i = 0;
                        }
                    }
                    matrices.add(m);
                }
            }
            for (JsonNode oneSeqAli : oneResult.get("sequence_alignment")) {
                seqs.add(oneSeqAli.get("sequence").asText());
//                for (JsonNode r : oneSeqAli.get("regions")) {
//                    System.out.println(r.toString());
//                }
            }
        }
        return new StrucAlignment(matrices, seqs);
    }

    private StrucAlignment computeFatcatAlignment(Chain c1, Chain c2) throws StructureException {
        // Get StructureAlignment instance
        StructureAlignment algorithm  = StructureAlignmentFactory.getAlgorithm(FatCatRigid.algorithmName);

        Atom[] ca1 = StructureTools.getAtomCAArray(c1);
        Atom[] ca2 = StructureTools.getAtomCAArray(c2);

        // Perform the alignment
        AFPChain afpChain = algorithm.align(ca1,ca2);

        // Print text output
        System.out.println(afpChain.toFatcat(ca1,ca2));

        // assuming 1 block
        Matrix m = afpChain.getBlockRotationMatrix()[0];
        Matrix4d transform = new Matrix4d();
        for (int i=0; i<3; i++) {
            for (int j=0; j<3; j++) {
                transform.setElement(i, j, m.get(j, i));
            }
        }
        Atom[] translation = afpChain.getBlockShiftVector();
        // assming 1 block
        transform.setTranslation(new Vector3d(translation[0].getX(), translation[0].getY(), translation[0].getZ()));
        transform.setElement(3,3,1);

        List<Matrix4d> transforms = Arrays.asList(ID_MATRIX, transform);
        List<String> seqs = Arrays.asList(new String(afpChain.getAlnseq1()), new String(afpChain.getAlnseq2()));

        return new StrucAlignment(transforms, seqs);
    }
}
